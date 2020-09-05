package seakers.engineerserver.server;

/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements. See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership. The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */


import java.io.File;
import java.util.*;
import java.util.concurrent.*;

import com.rabbitmq.client.Channel;
import com.rabbitmq.client.Connection;
import com.rabbitmq.client.ConnectionFactory;
import jess.Fact;
import jess.JessException;
import jess.Value;
import jess.ValueVector;
import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.real.UM;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.util.TypedProperties;

import seakers.architecture.util.IntegerVariable;
import seakers.engineerserver.search.BinaryInputInteractiveSearch;
import seakers.architecture.operators.IntegerUM;
import seakers.engineerserver.search.DiscreteInputInteractiveSearch;
import seakers.engineerserver.search.problems.Assigning.AssigningArchitecture;
import seakers.engineerserver.search.problems.Assigning.AssigningProblem;
import seakers.engineerserver.search.problems.PartitioningAndAssigning.PartitioningAndAssigningArchitecture;
import seakers.engineerserver.search.problems.PartitioningAndAssigning.PartitioningAndAssigningInitialization;
import seakers.engineerserver.thriftinterface.*;
import seakers.engineerserver.search.problems.PartitioningAndAssigning.PartitioningAndAssigningProblem;
import seakers.engineerserver.search.problems.PartitioningAndAssigning.operators.PartitioningAndAssigningCrossover;
import seakers.engineerserver.search.problems.PartitioningAndAssigning.operators.PartitioningAndAssigningMutation;
import seakers.engineerserver.search.problems.Scheduling.SchedulingProblem;
import seakers.engineerserver.search.problems.Scheduling.SchedulingArchitecture;
import seakers.engineerserver.search.problems.Scheduling.SchedulingInitialization;
import seakers.orekit.util.OrekitConfig;
import seakers.vassar.BaseParams;
import seakers.architecture.enumeration.*;
import seakers.architecture.enumeration.FullFactorial;
import seakers.architecture.pattern.*;
import seakers.architecture.pattern.Permuting;
import seakers.vassar.Resource;
import seakers.vassar.Result;
import seakers.vassar.architecture.AbstractArchitecture;
import seakers.vassar.evaluation.AbstractArchitectureEvaluator;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.problems.Assigning.*;
import seakers.vassar.problems.Scheduling.*;
import seakers.vassar.problems.PartitioningAndAssigning.Decadal2017AerosolsParams;

public class VASSARInterfaceHandler implements VASSARInterface.Iface {

    private String resourcesPath;
    private Map<String, BaseParams> paramsMap;
    private Map<String, ArchitectureEvaluationManager> architectureEvaluationManagerMap;
    private Map<String, Thread> gaThreads;

    public VASSARInterfaceHandler() {
        // Set a path to the project folder
        this.resourcesPath = "../VASSAR_resources";
        this.paramsMap = new HashMap<>();
        this.architectureEvaluationManagerMap = new HashMap<>();
        this.gaThreads = new HashMap<>();

        OrekitConfig.init(1, this.resourcesPath + File.separator + "orekit");
    }

    private Runnable generateBinaryInputGATask(String problem, List<BinaryInputArchitecture> dataset, String id) {
        return () -> {
            System.out.println("Starting GA for binary input data");

            ExecutorService pool = Executors.newFixedThreadPool(1);
            CompletionService<Algorithm> ecs = new ExecutorCompletionService<>(pool);

            //parameters and operators for seakers.vassar_server.search
            TypedProperties properties = new TypedProperties();
            //seakers.vassar_server.search paramaters set here
            int maxEvals = 3000;
            properties.setInt("maxEvaluations", maxEvals);
            properties.setInt("populationSize", dataset.size());
            double crossoverProbability = 1.0;
            properties.setDouble("crossoverProbability", crossoverProbability);
            double mutationProbability = 1. / 60.;
            properties.setDouble("mutationProbability", mutationProbability);
            Variation singlecross;
            Variation bitFlip;
            Variation intergerMutation;
            Initialization initialization;

            //setup for epsilon MOEA
            double[] epsilonDouble = new double[]{0.001, 1};

            //initialize problem
            BaseParams params = this.getProblemParameters(problem);
            ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);

            Problem assignmentProblem = new AssigningProblem(new int[]{1}, problem, AEM, params);

            // Create a solution for each input arch in the dataset
            List<Solution> initial = new ArrayList<>(dataset.size());
            for (BinaryInputArchitecture arch : dataset) {
                AssigningArchitecture new_arch = new AssigningArchitecture(new int[]{1},
                        params.getNumInstr(), params.getNumOrbits(), 2);

                for (int j = 1; j < new_arch.getNumberOfVariables(); ++j) {
                    BinaryVariable var = new BinaryVariable(1);
                    var.set(0, arch.inputs.get(j-1));
                    new_arch.setVariable(j, var);
                }
                new_arch.setObjective(0, -arch.outputs.get(0));
                new_arch.setObjective(1, arch.outputs.get(1));
                new_arch.setAlreadyEvaluated(true);
                initial.add(new_arch);
            }

            initialization = new InjectedInitialization(assignmentProblem, dataset.size(), initial);

            //initialize population structure for algorithm
            Population population = new Population();
            EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
            ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
            TournamentSelection selection = new TournamentSelection(2, comp);

            singlecross = new OnePointCrossover(crossoverProbability);
            bitFlip = new BitFlip(mutationProbability);
            intergerMutation = new IntegerUM(mutationProbability);
            CompoundVariation var = new CompoundVariation(singlecross, bitFlip, intergerMutation);

            Algorithm eMOEA = new EpsilonMOEA(assignmentProblem, population, archive, selection, var, initialization);
            ecs.submit(new BinaryInputInteractiveSearch(eMOEA, properties, id));

            // Message queue
            // Notify listeners of GA starting in username channel
            ConnectionFactory factory = new ConnectionFactory();
            factory.setHost(System.getenv("RABBITMQ_HOST"));
            String sendbackQueueName = id + "_gabrain";

            try (Connection connection = factory.newConnection(); Channel channel = connection.createChannel()) {
                channel.queueDeclare(sendbackQueueName, false, false, false, null);
                String message = "{ \"type\": \"ga_started\" }";
                channel.basicPublish("", sendbackQueueName, null, message.getBytes("UTF-8"));
            }
            catch (Exception e) {
                e.printStackTrace();
            }

            try {
                Algorithm alg = ecs.take().get();
            } catch (InterruptedException | ExecutionException ex) {
                ex.printStackTrace();
            }

            // Notify listeners of new architectures in username channel
            try (Connection connection = factory.newConnection(); Channel channel = connection.createChannel()) {
                channel.queueDeclare(sendbackQueueName, false, false, false, null);
                String message = "{ \"type\": \"ga_done\" }";
                channel.basicPublish("", sendbackQueueName, null, message.getBytes("UTF-8"));
            }
            catch (Exception e) {
                e.printStackTrace();
            }

            pool.shutdown();

            System.out.println("DONE");
        };
    }

    private Runnable generateDiscreteInputGATask(String problem, List<DiscreteInputArchitecture> dataset, String id) {
        return () -> {
            System.out.println("Starting GA for discrete input data");

            ExecutorService pool = Executors.newFixedThreadPool(1);
            CompletionService<Algorithm> ecs = new ExecutorCompletionService<>(pool);

            //parameters and operators for seakers.vassar_server.search
            TypedProperties properties = new TypedProperties();
            //search paramaters set here
            int maxEvals = 3000;
            properties.setInt("maxEvaluations", maxEvals);
            properties.setInt("populationSize", dataset.size());

            double crossoverProbability = 1.0;
            properties.setDouble("crossoverProbability", crossoverProbability);
            double mutationProbability = 1. / 6.;
            properties.setDouble("mutationProbability", mutationProbability);

            //setup for epsilon MOEA
            double[] epsilonDouble = new double[]{0.001, 1};

            //initialize problem
            BaseParams params = this.getProblemParameters(problem);
            ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);

            Problem partitioningAndAssigningProblem = new PartitioningAndAssigningProblem(problem, AEM, params);

            // Create a solution for each input arch in the dataset
            List<Solution> initial = new ArrayList<>(dataset.size());
            for (DiscreteInputArchitecture arch : dataset) {
                PartitioningAndAssigningArchitecture new_arch = new PartitioningAndAssigningArchitecture(
                        params.getNumInstr(), params.getNumOrbits(), 2);

                int numPartitioningVariables = params.getNumInstr();
                int numAssignmentVariables = params.getNumInstr();

                for (int j = 0; j < numPartitioningVariables; ++j) {
                    IntegerVariable var = new IntegerVariable(arch.inputs.get(j), 0, params.getNumInstr());
                    new_arch.setVariable(j, var);
                }

                for (int j = numPartitioningVariables; j < numPartitioningVariables + numAssignmentVariables; ++j) {
                    IntegerVariable var = new IntegerVariable(arch.inputs.get(j), -1, params.getNumOrbits());
                    new_arch.setVariable(j, var);
                }

                new_arch.setObjective(0, -arch.outputs.get(0));
                new_arch.setObjective(1, arch.outputs.get(1));
                new_arch.setAlreadyEvaluated(true);
                initial.add(new_arch);
            }

            Initialization initialization = new PartitioningAndAssigningInitialization(partitioningAndAssigningProblem,
                    dataset.size(), initial, params);

            //initialize population structure for algorithm
            Population population = new Population();
            EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
            ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
            TournamentSelection selection = new TournamentSelection(2, comp);

            // Define operators
            Variation singlecross = new PartitioningAndAssigningCrossover(crossoverProbability, params);
            Variation intergerMutation = new PartitioningAndAssigningMutation(mutationProbability, params);
            CompoundVariation var = new CompoundVariation(singlecross, intergerMutation);

            Algorithm eMOEA = new EpsilonMOEA(partitioningAndAssigningProblem, population, archive, selection, var, initialization);
            ecs.submit(new DiscreteInputInteractiveSearch(eMOEA, properties, id));

            // Message queue
            // Notify listeners of GA starting in username channel
            ConnectionFactory factory = new ConnectionFactory();
            factory.setHost(System.getenv("RABBITMQ_HOST"));
            String sendbackQueueName = id + "_gabrain";

            try (Connection connection = factory.newConnection(); Channel channel = connection.createChannel()) {
                channel.queueDeclare(sendbackQueueName, false, false, false, null);
                String message = "{ \"type\": \"ga_started\" }";
                channel.basicPublish("", sendbackQueueName, null, message.getBytes("UTF-8"));
            }
            catch (Exception e) {
                e.printStackTrace();
            }

            try {
                Algorithm alg = ecs.take().get();
            } catch (InterruptedException | ExecutionException ex) {
                ex.printStackTrace();
            }

            // Notify listeners of new architectures in username channel
            try (Connection connection = factory.newConnection(); Channel channel = connection.createChannel()) {
                channel.queueDeclare(sendbackQueueName, false, false, false, null);
                String message = " \"{ \"type\": \"ga_done\" }\"";
                channel.basicPublish("", sendbackQueueName, null, message.getBytes("UTF-8"));
            }
            catch (Exception e) {
                e.printStackTrace();
            }

            pool.shutdown();

            System.out.println("DONE");
        };
    }

    private Runnable generateSchedulingGATask(String problem, List<SchedulingInputArchitecture> dataset, List<BinaryInputArchitecture> inputArches, List<MissionMeasurements> historical_info, String id) {
        return () -> {
            System.out.println("Starting GA for scheduling");

            ExecutorService pool = Executors.newFixedThreadPool(1);
            CompletionService<Algorithm> ecs = new ExecutorCompletionService<>(pool);

            //parameters and operators for seakers.vassar_server.search
            TypedProperties properties = new TypedProperties();
            //search paramaters set here
            int maxEvals = 3000;
            int populationSize = 100;
            properties.setInt("maxEvaluations", maxEvals);
            properties.setInt("populationSize", populationSize);

            double crossoverProbability = 1.0;
            properties.setDouble("crossoverProbability", crossoverProbability);
            double mutationProbability = 1. / 6.;
            properties.setDouble("mutationProbability", mutationProbability);

            //setup for epsilon MOEA
            double[] epsilonDouble = new double[]{0.001, 0.001};

            //initialize problem
            BaseParams params = this.getProblemParameters(problem);

            Problem schedulingProblem = new SchedulingProblem(problem, this, inputArches.size(), inputArches, historical_info);

            // Create a solution for each input arch in the dataset
            List<Solution> initial = new ArrayList<>(dataset.size());
            for (SchedulingInputArchitecture arch : dataset) {
                SchedulingArchitecture new_arch = new SchedulingArchitecture(
                        params.getNumInstr(), 2);

                int numSchedulingVariables = inputArches.size();

                for (int j = 0; j < numSchedulingVariables; ++j) {
                    IntegerVariable var = new IntegerVariable(arch.inputs.get(j), 0, inputArches.size());
                    new_arch.setVariable(j, var);
                }

                new_arch.setObjective(0, -arch.outputs.get(0));
                new_arch.setObjective(1, -arch.outputs.get(1));
                new_arch.setAlreadyEvaluated(true);
                initial.add(new_arch);
            }

            Initialization initialization = new SchedulingInitialization(schedulingProblem,
                    dataset.size(), initial, params);

            //initialize population structure for algorithm
            Population population = new Population();
            EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
            ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
            TournamentSelection selection = new TournamentSelection(2, comp);

            // Define operators
            //Variation singlecross = new PartitioningAndAssigningCrossover(crossoverProbability, params);
            //Variation intergerMutation = new PartitioningAndAssigningMutation(mutationProbability, params);
            UM var = new UM(0.01);

            Algorithm eMOEA = new EpsilonMOEA(schedulingProblem, population, archive, selection, var, initialization);
            ecs.submit(new DiscreteInputInteractiveSearch(eMOEA, properties, id));

            // Message queue
            // Notify listeners of GA starting in username channel
            ConnectionFactory factory = new ConnectionFactory();
            factory.setHost(System.getenv("RABBITMQ_HOST"));
            String sendbackQueueName = id + "_gabrain";

            try (Connection connection = factory.newConnection(); Channel channel = connection.createChannel()) {
                channel.queueDeclare(sendbackQueueName, false, false, false, null);
                String message = "{ \"type\": \"ga_started\" }";
                channel.basicPublish("", sendbackQueueName, null, message.getBytes("UTF-8"));
            }
            catch (Exception e) {
                e.printStackTrace();
            }

            try {
                Algorithm alg = ecs.take().get();
            } catch (InterruptedException | ExecutionException ex) {
                ex.printStackTrace();
            }

            // Notify listeners of new architectures in username channel
            try (Connection connection = factory.newConnection(); Channel channel = connection.createChannel()) {
                channel.queueDeclare(sendbackQueueName, false, false, false, null);
                String message = " \"{ \"type\": \"ga_done\" }\"";
                channel.basicPublish("", sendbackQueueName, null, message.getBytes("UTF-8"));
            }
            catch (Exception e) {
                e.printStackTrace();
            }

            pool.shutdown();

            System.out.println("DONE");
        };
    }

    public void ping() {
        System.out.println("ping()");
    }

    private void initJess(String problem) {
        BaseParams params;
        AbstractArchitectureEvaluator evaluator;
        ArchitectureEvaluationManager AEM;
        String key;
        String search_clps = "";

        if (problem.equalsIgnoreCase("SMAP") || problem.equalsIgnoreCase("SMAP_JPL1")
                || problem.equalsIgnoreCase("SMAP_JPL2")
                || problem.equalsIgnoreCase("ClimateCentric")) {
            key = problem;

            if (problem.equalsIgnoreCase("SMAP")) {
                params = new SMAPParams(this.resourcesPath, "CRISP-ATTRIBUTES",
                        "test", "normal");
            }
            else if (problem.equalsIgnoreCase("SMAP_JPL1")) {
                params = new SMAPJPL1Params(this.resourcesPath, "CRISP-ATTRIBUTES",
                        "test", "normal");
            }
            else if (problem.equalsIgnoreCase("SMAP_JPL2")) {
                params = new SMAPJPL2Params(this.resourcesPath, "CRISP-ATTRIBUTES",
                        "test", "normal");
            }
            else if (problem.equalsIgnoreCase("ClimateCentric")) {
                params = new ClimateCentricParams(this.resourcesPath, "CRISP-ATTRIBUTES",
                        "test", "normal");
            }
            else {
                throw new RuntimeException();
            }

            evaluator = new seakers.vassar.problems.Assigning.ArchitectureEvaluator();

        }
        else if (problem.equalsIgnoreCase("Decadal2017Aerosols")) {
            key = "Decadal2017Aerosols";

            params = new Decadal2017AerosolsParams(this.resourcesPath, "CRISP-ATTRIBUTES",
                    "test", "normal");
            evaluator = new seakers.vassar.problems.PartitioningAndAssigning.ArchitectureEvaluator();
        }
        else {
            throw new IllegalArgumentException("Unrecorgnizable problem type: " + problem);
        }

        if (this.paramsMap.keySet().contains(key)) {
            // Already initialized
            return;
        }

        AEM = new ArchitectureEvaluationManager(params, evaluator);
        this.paramsMap.put(key, params);
        this.architectureEvaluationManagerMap.put(key, AEM);

        // Initialization
        AEM.init(2);
    }

    private BaseParams getProblemParameters(String problem) {
        String key;

        if (problem.equalsIgnoreCase("SMAP")) {
            key = "SMAP";
        }
        else if (problem.equalsIgnoreCase("SMAP_JPL1")) {
            key = "SMAP_JPL1";
        }
        else if (problem.equalsIgnoreCase("SMAP_JPL2")) {
            key = "SMAP_JPL2";
        }
        else if(problem.equalsIgnoreCase("ClimateCentric")){
            key = "ClimateCentric";
        }
        else if (problem.equalsIgnoreCase("Decadal2017Aerosols")) {
            key = "Decadal2017Aerosols";
        }
        else {
            throw new IllegalArgumentException("Unrecorgnizable problem type: " + problem);
        }

        if (this.paramsMap.containsKey(key)) {
            // Already initialized
            return this.paramsMap.get(problem);
        } else {
            this.initJess(problem);
            return this.paramsMap.get(problem);
        }
    }

    private AbstractArchitecture getArchitectureBinaryInput(String problem, String bitString, int numSatellites, BaseParams params) {
        AbstractArchitecture architecture;

        if (problem.equalsIgnoreCase("SMAP") || problem.equalsIgnoreCase("SMAP_JPL1")
                || problem.equalsIgnoreCase("SMAP_JPL2")
                || problem.equalsIgnoreCase("ClimateCentric")) {
            // Generate a new architecture
            architecture = new seakers.vassar.problems.Assigning.Architecture(bitString, numSatellites, (AssigningParams) params);

        } else {
            throw new IllegalArgumentException("Unrecorgnizable problem type: " + problem);
        }

        return architecture;
    }

    private AbstractArchitecture getArchitectureDiscreteInput(String problem, int[] intArray, int numSatellites, BaseParams params) {
        AbstractArchitecture architecture;

        if (problem.equalsIgnoreCase("Decadal2017Aerosols")) {
            // Generate a new architecture
            int numInstr = params.getNumInstr();
            int[] instrPartitioning = Arrays.copyOfRange(intArray, 0, numInstr);
            int[] orbitAssignation = Arrays.copyOfRange(intArray, numInstr, intArray.length);
            architecture = new seakers.vassar.problems.PartitioningAndAssigning.Architecture(instrPartitioning, orbitAssignation, numSatellites, (Decadal2017AerosolsParams) params);

        } else {
            throw new IllegalArgumentException("Unrecorgnizable problem type: " + problem);
        }

        return architecture;
    }

    @Override
    public double evaluateDataContinuityScore(List<MissionMeasurements> missionMeasurements, List<MissionMeasurements> historical_missionMeasurements) {
        int score = 0;
        for(MissionMeasurements missionMeasurement : missionMeasurements) {
            int start_year = missionMeasurement.start_year;
            int end_year = missionMeasurement.end_year;
            for(String individualMeasurement : missionMeasurement.names) {
                for(MissionMeasurements historical_missionMeasurement : historical_missionMeasurements) {
                    for(String name : historical_missionMeasurement.names) {
                        if(name.equals("Sea-ice cover")) {
                            name = "Sea ice cover";
                        }
                        if(name.equals("Soil moisture at the surface")) {
                            name = "Soil moisture";
                        }
                        if(name.equals("Wind vector over sea surface (horizontal)") || name.equals("Wind speed over sea surface (horizontal)")) {
                            name = "Ocean surface wind speed";
                        }
                        if(name.equals(individualMeasurement)) {
                            if (end_year < historical_missionMeasurement.start_year || start_year > historical_missionMeasurement.end_year) {
                                score = score;
                            } else if (start_year > historical_missionMeasurement.start_year && end_year < historical_missionMeasurement.end_year) {
                                score = score + end_year - start_year;
                            } else if (start_year > historical_missionMeasurement.start_year && end_year > historical_missionMeasurement.end_year) {
                                score = score + historical_missionMeasurement.end_year - start_year;
                            } else if (start_year < historical_missionMeasurement.start_year && end_year < historical_missionMeasurement.end_year) {
                                score = score + end_year - historical_missionMeasurement.start_year;
                            } else {
                                score = score;
                            }
                            //System.out.println("Score: "+score+", start_year: "+start_year+", end_year: "+ end_year+", missionMeasurement.start_year: "+missionMeasurement.start_year+", missionMeasurement.end_year: "+missionMeasurement.end_year);
                        }
                    }
                }
            }
        }
        return Double.valueOf(score);
    }

    @Override
    public double evaluateFairnessScore(List<MissionMeasurements> missionMeasurements) {
        int score = 0;
        List<Double> panelSum = new ArrayList<Double>(Arrays.asList(0.0,0.0,0.0,0.0,0.0));
        for (MissionMeasurements missionMeasurement : missionMeasurements) {
            List<Double> panelScores = missionMeasurement.panelScores;
            for (int i = 0; i < panelScores.size(); i++) {
                panelSum.set(i,panelSum.get(i) + panelScores.get(i));
            }
        }
        double minimum = 100000.0;
        for (int j = 0; j < panelSum.size(); j++) {
            if(panelSum.get(j)<minimum) {
                minimum = panelSum.get(j);
            }
        }
        return minimum;
    }

    @Override
    public double enumeratedDesigns(String problem, List<BinaryInputArchitecture> input_arches, List<MissionMeasurements> historical_info) {
        double xd = 250.0;
        double lowestScore = 100000.0;
        int startYear = 2010;
        int yearIncrement = 2;
        int numArches = input_arches.size();
        int[] lowestOption = new int[numArches];
        Permuting permute = new Permuting(numArches,"xd");
        FullFactorial ff = new FullFactorial(permute);
        Collection<int[]> options = ff.ffPermuting(numArches);
        for (int[] option : options) {
            //System.out.println(Arrays.toString(option));
            List<MissionMeasurements> missions = new ArrayList<>();
            for (int i = 0; i < option.length; i++) {
                int missionStartYear = startYear + option[i]*yearIncrement;
                int missionEndYear = missionStartYear + yearIncrement;
                //System.out.println("Start Year: "+missionStartYear+", End Year: "+missionEndYear);
                List<String> measurements = getMeasurements(problem, input_arches.get(i));
                List<Double> panelScores = getPanelScoresForArch(problem, input_arches.get(i));
                //System.out.println(Arrays.toString(measurements.toArray()));
                MissionMeasurements mission = new MissionMeasurements(measurements, missionStartYear, missionEndYear, panelScores);
                missions.add(mission);
            }
            double dataContinuityScore = evaluateDataContinuityScore(missions, historical_info);
            if(dataContinuityScore < lowestScore) {
                lowestScore = dataContinuityScore;
                lowestOption = option;
            }
            double fairnessScore = evaluateFairnessScore(missions);
            System.out.println("Data Continuity Score: " + dataContinuityScore + ", Fairness Score: " + fairnessScore);
        }
        System.out.println("Best order: "+Arrays.toString(lowestOption));
        return xd;
    }

    @Override
    public BinaryInputArchitecture evalBinaryInputArch(String problem, List<Boolean> boolList) {
        // Input a new architecture design
        String bitString = "";
        for (Boolean b : boolList) {
            bitString += b ? "1" : "0";
        }

        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        AbstractArchitecture architecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);

        // Evaluate the architecture
        Result result = AEM.evaluateArchitectureSync(architecture, "Slow");

        // Save the score and the cost
        double cost = result.getCost();
        double science = result.getScience();
        List<Double> outputs = new ArrayList<>();
        outputs.add(science);
        outputs.add(cost);

        System.out.println("Performance Score: " + science + ", Cost: " + cost);
        return new BinaryInputArchitecture(0, boolList, outputs);
    }

    @Override
    public DiscreteInputArchitecture evalDiscreteInputArch(String problem, List<Integer> intList) {
        // Input a new architecture design
        int[] intArray = new int[intList.size()];
        for (int i = 0; i < intList.size(); i++) {
            intArray[i] = intList.get(i);
        }

        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        AbstractArchitecture architecture = this.getArchitectureDiscreteInput(problem, intArray, 1, params);

        // Evaluate the architecture
        Result result = AEM.evaluateArchitectureSync(architecture, "Slow");

        // Save the score and the cost
        double cost = result.getCost();
        double science = result.getScience();
        List<Double> outputs = new ArrayList<>();
        outputs.add(science);
        outputs.add(cost);

        System.out.println("Performance Score: " + science + ", Cost: " + cost);
        return new DiscreteInputArchitecture(0, intList, outputs);
    }


    @Override
    public List<BinaryInputArchitecture> runLocalSearchBinaryInput(String problem, BinaryInputArchitecture inputArch) {

        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);

        String bitString = "";
        for (Boolean b : inputArch.inputs) {
            bitString += b ? "1" : "0";
        }

        ArrayList<String> samples = randomLocalChangeBinaryInput(bitString, 4, params);

        List<BinaryInputArchitecture> out = new ArrayList<>();
        for (String sample : samples) {

            AbstractArchitecture architecture = this.getArchitectureBinaryInput(problem, sample, 1, params);

            // Evaluate the architecture
            Result result = AEM.evaluateArchitectureSync(architecture, "Slow");

            // Save the score and the cost
            double cost = result.getCost();
            double science = result.getScience();
            List<Double> outputs = new ArrayList<>();
            outputs.add(science);
            outputs.add(cost);

            System.out.println("bitString: " + sample + ", Science: " + science + ", Cost: " + cost);

            BinaryInputArchitecture arch = new BinaryInputArchitecture(0, bitString2BoolArray(sample), outputs);
            out.add(arch);
        }

        return out;
    }

    @Override
    public List<DiscreteInputArchitecture> runLocalSearchDiscreteInput(String problem, DiscreteInputArchitecture inputArch) {
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);

        int[] inputsArray = new int[inputArch.inputs.size()];
        for(int i = 0; i < inputArch.inputs.size(); i++){
            inputsArray[i] = inputArch.inputs.get(i);
        }

        ArrayList<int[]> samples = randomLocalChangeDiscreteInput(inputsArray, 4, params);

        List<DiscreteInputArchitecture> out = new ArrayList<>();
        for (int[] sample : samples) {

            AbstractArchitecture architecture = this.getArchitectureDiscreteInput(problem, sample, 1, params);

            // Evaluate the architecture
            Result result = AEM.evaluateArchitectureSync(architecture, "Slow");

            // Save the score and the cost
            double cost = result.getCost();
            double science = result.getScience();
            List<Double> outputs = new ArrayList<>();
            outputs.add(science);
            outputs.add(cost);

            List<Integer> newInput = new ArrayList<>();
            for(int i = 0; i < sample.length; i++){
                newInput.add(sample[i]);
            }

            DiscreteInputArchitecture arch = new DiscreteInputArchitecture(0, newInput, outputs);
            out.add(arch);
        }

        return out;
    }

    private ArrayList<int[]> randomLocalChangeDiscreteInput(int[] input, int n, BaseParams params) {

        // Assumes partitioning + assigning problem

        Random rand = new Random();
        int numVars = 2 * params.getNumInstr();

        ArrayList<int[]> out = new ArrayList<>();
        Solution solution = new PartitioningAndAssigningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2);

        for(int i = 0; i < n; i++){
            int selectedVar = rand.nextInt(numVars);
            IntegerVariable variable = (IntegerVariable) solution.getVariable(selectedVar);
            int lowerBound = variable.getLowerBound();
            int upperBound = variable.getUpperBound();

            int randInt = rand.nextInt(upperBound - lowerBound + 1) + lowerBound;
            int[] modifiedInput = Arrays.copyOfRange(input, 0, input.length);
            modifiedInput[selectedVar] = randInt;

            out.add(modifiedInput);
        }

        return out;
    }

    private ArrayList<String> randomLocalChangeBinaryInput(String bitString, int n, BaseParams params) {
        Random rand = new Random();
        int numVars = params.getNumOrbits() * params.getNumInstr();

        ArrayList<String> out = new ArrayList<>();

        for (int i = 0; i < n; i++) {
            int k = rand.nextInt(numVars);

            StringBuilder tempBitString = new StringBuilder(bitString);
            if (bitString.charAt(k) == '1') {
                tempBitString.setCharAt(k, '0');
            } else {
                tempBitString.setCharAt(k, '1');
            }
            out.add(tempBitString.toString());
        }
        return out;
    }

    private List<Boolean> bitString2BoolArray(String bitString) {
        List<Boolean> out = new ArrayList<>();
        for (int i = 0; i < bitString.length(); i++) {
            out.add(bitString.charAt(i) == '1');
        }
        return out;
    }

    @Override
    public List<String> getCritiqueBinaryInputArch(String problem, BinaryInputArchitecture inputArch) {

        String bitString = "";
        for (Boolean b : inputArch.inputs) {
            bitString += b ? "1" : "0";
        }

        AbstractArchitecture architecture;

        System.out.println(bitString);

        if (problem.equalsIgnoreCase("SMAP") || problem.equalsIgnoreCase("SMAP_JPL1")
                || problem.equalsIgnoreCase("SMAP_JPL2")
                || problem.equalsIgnoreCase("ClimateCentric")) {

            seakers.vassar.problems.Assigning.AssigningParams params = (seakers.vassar.problems.Assigning.AssigningParams) this.getProblemParameters(problem);
            ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);

            // Generate a new architecture
            architecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);

            // Initialize Critique Generator
            seakers.vassar.problems.Assigning.CritiqueGenerator critiquer = new seakers.vassar.problems.Assigning.CritiqueGenerator(AEM.getResourcePool(), architecture);

            return critiquer.getCritique();

        } else {
            throw new IllegalArgumentException("Unrecorgnizable problem type: " + problem);
        }
    }

    @Override
    public List<String> getCritiqueDiscreteInputArch(String problem, DiscreteInputArchitecture inputArch) {

        AbstractArchitecture architecture;

        System.out.println(inputArch.inputs);

        if (problem.equalsIgnoreCase("Decadal2017Aerosols")) {

            seakers.vassar.problems.PartitioningAndAssigning.PartitioningAndAssigningParams params = (seakers.vassar.problems.PartitioningAndAssigning.PartitioningAndAssigningParams) this.getProblemParameters(problem);
            ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);

            // Generate a new architecture
            int[] intArray = new int[inputArch.inputs.size()];
            for(int i = 0; i < intArray.length; i++)
                intArray[i] = inputArch.inputs.get(i);
            architecture = this.getArchitectureDiscreteInput(problem, intArray, 1, params);

            // Initialize Critique Generator
            seakers.vassar.problems.PartitioningAndAssigning.CritiqueGenerator critiquer = new seakers.vassar.problems.PartitioningAndAssigning.CritiqueGenerator(AEM.getResourcePool(), architecture);

            return critiquer.getCritique();

        } else {
            throw new IllegalArgumentException("Unrecorgnizable problem type: " + problem);
        }
    }

    @Override
    public ArrayList<String> getOrbitList(String problem) {
        this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        ArrayList<String> orbitList = new ArrayList<>();
        for (String o : params.getOrbitList()) {
            orbitList.add(o);
        }
        return orbitList;
    }

    @Override
    public ArrayList<String> getInstrumentList(String problem) {
        this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        ArrayList<String> instrumentList = new ArrayList<>();
        for (String i : params.getInstrumentList()) {
            instrumentList.add(i);
        }
        return instrumentList;
    }

    @Override
    public ArrayList<String> getObjectiveList(String problem) {
        this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        ArrayList<String> objectiveList = new ArrayList<>();
        params.objectiveDescriptions.forEach((k, v) -> {
            objectiveList.add(k);
        });
        return objectiveList;
    }

    @Override
    public ArrayList<String> getSubobjectiveList(String problem) {
        this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        ArrayList<String> subobjectiveList = new ArrayList<>();
        params.subobjDescriptions.forEach((k, v) -> {
            subobjectiveList.add(k);
        });
        return subobjectiveList;
    }

    @Override
    public ArrayList<String> getInstrumentsForObjective(String problem, String objective) {
        this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        return new ArrayList<>(params.objectivesToInstruments.get(objective));
    }

    @Override
    public ArrayList<String> getInstrumentsForPanel(String problem, String panel) {
        this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        return new ArrayList<>(params.panelsToInstruments.get(panel));
    }

    @Override
    public List<ObjectiveSatisfaction> getArchitectureScoreExplanation(String problem, BinaryInputArchitecture arch) {
        String bitString = "";
        for (Boolean b : arch.inputs) {
            bitString += b ? "1" : "0";
        }

        this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);

        // Generate a new architecture
        AbstractArchitecture architecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);

        // Evaluate the architecture
        Result result = null;
        // Save the explanations for each stakeholder score
        List<ObjectiveSatisfaction> explanations = new ArrayList<>();

        result = AEM.evaluateArchitectureSync(architecture, "Slow");
        for (int i = 0; i < params.panelNames.size(); ++i) {
            explanations.add(new ObjectiveSatisfaction(params.panelNames.get(i),
                    result.getPanelScores().get(i), params.panelWeights.get(i)));
        }

        return explanations;
    }


    @Override
    public List<ObjectiveSatisfaction> getPanelScoreExplanation(String problem, BinaryInputArchitecture arch, String panel) {
        String bitString = "";
        for (Boolean b : arch.inputs) {
            bitString += b ? "1" : "0";
        }

        this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);

        // Generate a new architecture
        AbstractArchitecture architecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);

        // Evaluate the architecture
        Result result = null;
        // Save the explanations for each stakeholder score
        List<ObjectiveSatisfaction> explanations = new ArrayList<>();

        result = AEM.evaluateArchitectureSync(architecture, "Slow");
        for (int i = 0; i < params.panelNames.size(); ++i) {
            if (params.panelNames.get(i).equals(panel)) {
                for (int j = 0; j < params.objNames.get(i).size(); ++j) {
                    explanations.add(new ObjectiveSatisfaction(params.objNames.get(i).get(j),
                            result.getObjectiveScores().get(i).get(j), params.objWeights.get(i).get(j)));
                }
            }
        }

        return explanations;
    }

    @Override
    public List<ObjectiveSatisfaction> getObjectiveScoreExplanation(String problem, BinaryInputArchitecture arch, String objective) {
        String bitString = "";
        for (Boolean b : arch.inputs) {
            bitString += b ? "1" : "0";
        }

        this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);

        // Generate a new architecture
        AbstractArchitecture architecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);

        // Evaluate the architecture
        Result result = null;
        // Save the explanations for each stakeholder score
        List<ObjectiveSatisfaction> explanations = new ArrayList<>();

        result = AEM.evaluateArchitectureSync(architecture, "Slow");
        for (int i = 0; i < params.panelNames.size(); ++i) {
            for (int j = 0; j < params.objNames.get(i).size(); ++j) {
                if (params.objNames.get(i).get(j).equals(objective)) {
                    for (int k = 0; k < params.subobjectives.get(i).get(j).size(); ++k) {
                        explanations.add(new ObjectiveSatisfaction(params.subobjectives.get(i).get(j).get(k),
                                result.getSubobjectiveScores().get(i).get(j).get(k),
                                params.subobjWeights.get(i).get(j).get(k)));
                    }
                }
            }
        }

        return explanations;
    }

    private SubobjectiveDetails getSubscoreDetails(BaseParams params, String subobj, Result result) {
        String parameter = params.subobjectivesToMeasurements.get(subobj);

        // Obtain list of attributes for this parameter
        ArrayList<String> attrNames = new ArrayList<>();
        HashMap<String, ArrayList<String>> requirementRules = params.requirementRules.get(subobj);
        attrNames.addAll(requirementRules.keySet());
        HashMap<String, Integer> numDecimals = new HashMap<>();
        numDecimals.put("Horizontal-Spatial-Resolution#", 0);
        numDecimals.put("Temporal-resolution#", 0);
        numDecimals.put("Swath#", 0);

        // Loop to get rows of details for each data product
        ArrayList<List<String>> attrValues = new ArrayList<>();
        ArrayList<Double> scores = new ArrayList<>();
        ArrayList<String> takenBy = new ArrayList<>();
        ArrayList<List<String>> justifications = new ArrayList<>();
        for (Fact explanation: result.getExplanations().get(subobj)) {
            try {
                // Try to find the requirement fact!
                int measurementId = explanation.getSlotValue("requirement-id").intValue(null);
                if (measurementId == -1) {
                    continue;
                }
                Fact measurement = null;
                for (Fact capability: result.getCapabilities()) {
                    if (capability.getFactId() == measurementId) {
                        measurement = capability;
                        break;
                    }
                }
                // Start by putting all attribute values into list
                ArrayList<String> rowValues = new ArrayList<>();
                for (String attrName: attrNames) {
                    String attrType = requirementRules.get(attrName).get(0);
                    // Check type and convert to String if needed
                    Value attrValue = measurement.getSlotValue(attrName);
                    switch (attrType) {
                        case "SIB":
                        case "LIB": {
                            Double value = attrValue.floatValue(null);
                            double scale = 100;
                            if (numDecimals.containsKey(attrName)) {
                                scale = Math.pow(10, numDecimals.get(attrName));
                            }
                            value = Math.round(value * scale) / scale;
                            rowValues.add(value.toString());
                            break;
                        }
                        default: {
                            rowValues.add(attrValue.toString());
                            break;
                        }
                    }
                }
                // Get information from explanation fact
                Double score = explanation.getSlotValue("satisfaction").floatValue(null);
                String satisfiedBy = explanation.getSlotValue("satisfied-by").stringValue(null);
                ArrayList<String> rowJustifications = new ArrayList<>();
                ValueVector reasons = explanation.getSlotValue("reasons").listValue(null);
                for (int i = 0; i < reasons.size(); ++i) {
                    String reason = reasons.get(i).stringValue(null);
                    if (!reason.equals("N-A")) {
                        rowJustifications.add(reason);
                    }
                }

                // Put everything in their lists
                attrValues.add(rowValues);
                scores.add(score);
                takenBy.add(satisfiedBy);
                justifications.add(rowJustifications);
            }
            catch (JessException e) {
                System.err.println(e.toString());
            }
        }

        return new SubobjectiveDetails(
                parameter,
                attrNames,
                attrValues,
                scores,
                takenBy,
                justifications);
    }

    @Override
    public SubobjectiveDetails getSubscoreDetailsBinaryInput(String problem, BinaryInputArchitecture architecture, String subobj) {
        // Get a result with all the important facts
        String bitString = "";
        for (Boolean b : architecture.inputs) {
            bitString += b ? "1" : "0";
        }

        this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        AbstractArchitecture absArchitecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);
        Result result = AEM.evaluateArchitectureSync(absArchitecture, "Slow", true);

        return getSubscoreDetails(params, subobj, result);
    }

    @Override
    public SubobjectiveDetails getSubscoreDetailsDiscreteInput(String problem, DiscreteInputArchitecture architecture, String subobj) {
        // Get a result with all the important facts
        int[] intArray = new int[architecture.inputs.size()];
        for (int i = 0; i < architecture.inputs.size(); i++) {
            intArray[i] = architecture.inputs.get(i);
        }

        BaseParams params = this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        AbstractArchitecture absArchitecture = this.getArchitectureDiscreteInput(problem, intArray, 1, params);

        // Evaluate the architecture
        Result result = AEM.evaluateArchitectureSync(absArchitecture, "Slow", true);

        return getSubscoreDetails(params, subobj, result);
    }

    private List<SubscoreInformation> getArchScienceInformation(BaseParams params, Result result) {
        List<SubscoreInformation> information = new ArrayList<>();

        for (int i = 0; i < params.panelNames.size(); ++i) {
            List<SubscoreInformation> objectivesInformation = new ArrayList<>();
            for (int j = 0; j < params.objNames.get(i).size(); ++j) {
                List<SubscoreInformation> subobjectivesInformation = new ArrayList<>();
                for (int k = 0; k < params.subobjectives.get(i).get(j).size(); ++k) {
                    String subobjName = params.subobjectives.get(i).get(j).get(k);
                    subobjectivesInformation.add(new SubscoreInformation(
                            subobjName,
                            params.subobjDescriptions.get(subobjName),
                            result.getSubobjectiveScores().get(i).get(j).get(k),
                            params.subobjWeights.get(i).get(j).get(k),
                            null));
                }
                String objName = params.objNames.get(i).get(j);
                objectivesInformation.add(new SubscoreInformation(
                        objName,
                        params.objectiveDescriptions.get(objName),
                        result.getObjectiveScores().get(i).get(j),
                        params.objWeights.get(i).get(j),
                        subobjectivesInformation));
            }
            String panelName = params.panelNames.get(i);
            information.add(new SubscoreInformation(
                    panelName,
                    params.panelDescriptions.get(panelName),
                    result.getPanelScores().get(i),
                    params.panelWeights.get(i),
                    objectivesInformation));
        }

        return information;
    }

    @Override
    public List<String> getCommonMeasurements(String problem, List<BinaryInputArchitecture> arch_list) {
        List<String> commonMeasurements = new ArrayList<>();
        for (BinaryInputArchitecture arch : arch_list) {
            //System.out.println("New architecture");
            List<SubscoreInformation> measurements = getArchScienceInformationBinaryInput(problem, arch);
            for (SubscoreInformation measurement : measurements) {
                List<SubscoreInformation> panel_subscores = measurement.subscores;
                for(SubscoreInformation panel_subscore : panel_subscores) {
                    List<SubscoreInformation> obj_subscores = panel_subscore.subscores;
                    for(SubscoreInformation obj_subscore : obj_subscores) {
                        String name = obj_subscore.name;
                        String description = obj_subscore.description;
                        double value = obj_subscore.value;
                        if(value != 0 && !commonMeasurements.contains(name+" "+description)) {
                            commonMeasurements.add(description);
                            //System.out.println(name+" "+description);
                        }
                    }
                }
            }
        }
        //System.out.println(commonMeasurements);
        return commonMeasurements;
    }

    @Override
    public List<String> getMeasurements(String problem, BinaryInputArchitecture arch) {
        List<String> archMeasurements = new ArrayList<>();
        List<SubscoreInformation> measurements = getArchScienceInformationBinaryInput(problem, arch);
        for (SubscoreInformation measurement : measurements) {
            List<SubscoreInformation> panel_subscores = measurement.subscores;
            for(SubscoreInformation panel_subscore : panel_subscores) {
                List<SubscoreInformation> obj_subscores = panel_subscore.subscores;
                for(SubscoreInformation obj_subscore : obj_subscores) {
                    String description = obj_subscore.description;
                    if(obj_subscore.value != 0) {
                        archMeasurements.add(description);
                    }
                }
            }
        }
        return archMeasurements;
    }

    @Override
    public List<Double> getPanelScoresForArch(String problem, BinaryInputArchitecture arch) {
        List<Double> panelScores = new ArrayList<>();
        List<SubscoreInformation> measurements = getArchScienceInformationBinaryInput(problem, arch);
        for (SubscoreInformation measurement : measurements) {
            double score = measurement.value;
            panelScores.add(score);
        }
        return panelScores;
    }

    @Override
    public List<SubscoreInformation> getArchScienceInformationBinaryInput(String problem, BinaryInputArchitecture architecture) {
        String bitString = "";
        for (Boolean b : architecture.inputs) {
            bitString += b ? "1" : "0";
        }

        this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        AbstractArchitecture absArchitecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);
        Result result = AEM.evaluateArchitectureSync(absArchitecture, "Slow", true);

        return getArchScienceInformation(params, result);
    }

    @Override
    public List<SubscoreInformation> getArchScienceInformationDiscreteInput(String problem, DiscreteInputArchitecture architecture) {
        int[] intArray = new int[architecture.inputs.size()];
        for (int i = 0; i < architecture.inputs.size(); i++) {
            intArray[i] = architecture.inputs.get(i);
        }

        this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        AbstractArchitecture absArchitecture = this.getArchitectureDiscreteInput(problem, intArray, 1, params);

        // Evaluate the architecture
        Result result = AEM.evaluateArchitectureSync(absArchitecture, "Slow", true);

        return getArchScienceInformation(params, result);
    }

    private List<MissionCostInformation> getArchCostInformation(Result result) {
        List<MissionCostInformation> information = new ArrayList<>();

        // Auxiliary arrays
        String[] massBudgetSlots = { "adapter-mass", "propulsion-mass#", "structure-mass#", "avionics-mass#",
                "ADCS-mass#", "EPS-mass#", "propellant-mass-injection", "propellant-mass-ADCS", "thermal-mass#",
                "payload-mass#" };
        String[] powerBudgetSlots = { "payload-peak-power#", "satellite-BOL-power#" };
        String[] costBudgetSlots = { "payload-cost#", "bus-cost#", "launch-cost#", "program-cost#",
                "IAT-cost#", "operations-cost#" };
        double[] costMultipliers = { 1e-3, 1e-3, 1.0, 1e-3, 1e-3, 1e-3 };
        for (Fact costFact: result.getCostFacts()) {
            try {
                String missionName = costFact.getSlotValue("Name").stringValue(null);
                // Obtain the list of instruments for this orbit
                ArrayList<String> payloads = new ArrayList<>();
                ValueVector instruments = costFact.getSlotValue("instruments").listValue(null);
                for (int i = 0; i < instruments.size(); ++i) {
                    payloads.add(instruments.get(i).stringValue(null));
                }

                // Get the launch vehicle name
                String launchVehicle = costFact.getSlotValue("launch-vehicle").stringValue(null);
                HashMap<String, Double> massBudget = new HashMap<>();
                for (String massSlot: massBudgetSlots) {
                    Double value = costFact.getSlotValue(massSlot).floatValue(null);
                    massBudget.put(massSlot, value);
                }
                HashMap<String, Double> powerBudget = new HashMap<>();
                Double totalPower = 0.0;
                for (String powerSlot: powerBudgetSlots) {
                    Double value = costFact.getSlotValue(powerSlot).floatValue(null);
                    totalPower += value;
                    powerBudget.put(powerSlot, value);
                }
                HashMap<String, Double> costBudget = new HashMap<>();
                Double sumCost = 0.0;
                for (int i = 0; i < costBudgetSlots.length; ++i) {
                    String costSlot = costBudgetSlots[i];
                    Double multiplier = costMultipliers[i];
                    Double value = costFact.getSlotValue(costSlot).floatValue(null);
                    sumCost += value*multiplier;
                    costBudget.put(costSlot, value*multiplier);
                }
                Double totalCost = costFact.getSlotValue("mission-cost#").floatValue(null);
                costBudget.put("others", totalCost - sumCost);
                Double totalMass = costFact.getSlotValue("satellite-launch-mass").floatValue(null);
                information.add(new MissionCostInformation(
                        missionName,
                        payloads,
                        launchVehicle,
                        totalMass,
                        totalPower,
                        totalCost,
                        massBudget,
                        powerBudget,
                        costBudget));
            }
            catch (JessException e) {
                System.err.println(e.toString());
            }
        }

        return information;
    }

    @Override
    public List<MissionCostInformation> getArchCostInformationBinaryInput(String problem, BinaryInputArchitecture architecture) {
        String bitString = "";
        for (Boolean b : architecture.inputs) {
            bitString += b ? "1" : "0";
        }

        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        AbstractArchitecture absArchitecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);
        Result result = AEM.evaluateArchitectureSync(absArchitecture, "Slow", true);

        return getArchCostInformation(result);
    }

    @Override
    public List<MissionCostInformation> getArchCostInformationDiscreteInput(String problem, DiscreteInputArchitecture architecture) {
        int[] intArray = new int[architecture.inputs.size()];
        for (int i = 0; i < architecture.inputs.size(); i++) {
            intArray[i] = architecture.inputs.get(i);
        }

        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        Resource res = AEM.getResourcePool().getResource();
        BaseParams params = res.getParams();
        AEM.getResourcePool().freeResource(res);
        AbstractArchitecture absArchitecture = this.getArchitectureDiscreteInput(problem, intArray, 1, params);

        // Evaluate the architecture
        Result result = AEM.evaluateArchitectureSync(absArchitecture, "Slow", true);

        return getArchCostInformation(result);
    }

    private void deleteAllDoneGAs() {
        this.gaThreads.entrySet().removeIf(entry -> !entry.getValue().isAlive());
    }

    @Override
    public boolean isGARunning(String id) {
        return this.gaThreads.get(id).isAlive();
    }

    @Override
    public int stopGA(String id) {
        if (this.gaThreads.containsKey(id) && this.gaThreads.get(id).isAlive())  {
            ConnectionFactory factory = new ConnectionFactory();
            factory.setHost(System.getenv("RABBITMQ_HOST"));
            String queueName = id + "_brainga";

            try (Connection connection = factory.newConnection(); Channel channel = connection.createChannel()) {
                channel.queueDeclare(queueName, false, false, false, null);
                String message = "close";
                channel.basicPublish("", queueName, null, message.getBytes("UTF-8"));
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            return 0;
        }

        // Remove all dead Threads
        deleteAllDoneGAs();
        return 1;
    }

    @Override
    public String startGABinaryInput(String problem, List<BinaryInputArchitecture> dataset, String username) {
        String gaId = UUID.randomUUID().toString() + "_" + username + "_" + problem;
        deleteAllDoneGAs();
        this.gaThreads.put(gaId, new Thread(this.generateBinaryInputGATask(problem, dataset, gaId)));
        this.gaThreads.get(gaId).start();
        return gaId;
    }

    @Override
    public String startGADiscreteInput(String problem, List<DiscreteInputArchitecture> dataset, String username) {
        String gaId = UUID.randomUUID().toString() + "_" + username + "_" + problem;
        deleteAllDoneGAs();
        this.gaThreads.put(gaId, new Thread(this.generateDiscreteInputGATask(problem, dataset, username)));
        this.gaThreads.get(gaId).start();
        return gaId;
    }

    @Override
    public String startGAScheduling(String problem, List<SchedulingInputArchitecture> dataset, List<BinaryInputArchitecture> inputArches, List<MissionMeasurements> historicalInfo, String username) {
        String gaId = UUID.randomUUID().toString() + "_" + username + "_" + problem;
        deleteAllDoneGAs();
        this.gaThreads.put(gaId, new Thread(this.generateSchedulingGATask(problem, dataset, inputArches, historicalInfo, gaId)));
        this.gaThreads.get(gaId).start();
        return gaId;
    }
}
