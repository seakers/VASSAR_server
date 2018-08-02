package server;

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

import io.lettuce.core.RedisClient;
import javaInterface.*;
import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.util.TypedProperties;

import seak.architecture.util.IntegerVariable;
import search.InteractiveSearch;
import rbsa.eoss.architecture.AbstractArchitecture;
import rbsa.eoss.evaluation.AbstractArchitectureEvaluator;
import rbsa.eoss.evaluation.ArchitectureEvaluationManager;
import seak.architecture.operators.IntegerUM;
import rbsa.eoss.local.BaseParams;
import rbsa.eoss.Result;
import search.problems.Assigning.AssigningArchitecture;
import search.problems.Assigning.AssigningProblem;
import search.problems.PartitioningAndAssigning.PartitioningAndAssigningArchitecture;

public class VASSARInterfaceHandler implements VASSARInterface.Iface {

    private String root;
    private Map<String, BaseParams> paramsMap;
    private Map<String, ArchitectureEvaluationManager> architectureEvaluationManagerMap;

    public VASSARInterfaceHandler() {
        // Set a path to the project folder
        this.root = System.getProperty("user.dir");
        this.paramsMap = new HashMap<>();
        this.architectureEvaluationManagerMap = new HashMap<>();
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

        if (problem.equalsIgnoreCase("SMAP") ||
                problem.equalsIgnoreCase("ClimateCentric")) {

            key = problem;
            String path = this.root +
                    File.separator + "problems" +
                    File.separator + key;

            if(problem.equalsIgnoreCase("SMAP")){
                params = new rbsa.eoss.problems.Assigning.SMAPParams(path, "FUZZY-ATTRIBUTES", "test", "normal", search_clps);

            }else if(problem.equalsIgnoreCase("ClimateCentric")){
                params = new rbsa.eoss.problems.Assigning.ClimateCentricParams(path, "FUZZY-ATTRIBUTES", "test", "normal", search_clps);

            }else{
                throw new RuntimeException();
            }

            evaluator = new rbsa.eoss.problems.Assigning.ArchitectureEvaluator(params);

        } else if (problem.equalsIgnoreCase("Decadal2017Aerosols")) {

            key = "Decadal2017Aerosols";
            String path = this.root +
                    File.separator + "problems" +
                    File.separator + "SMAP";
            params = new rbsa.eoss.problems.PartitioningAndAssigning.Decadal2017AerosolsParams(path, "FUZZY-ATTRIBUTES", "test", "normal", search_clps);
            evaluator = new rbsa.eoss.problems.PartitioningAndAssigning.ArchitectureEvaluator(params);


        } else {
            throw new IllegalArgumentException("Unrecorgnizable problem type: " + problem);
        }

        if (this.paramsMap.keySet().contains(key)) {
            // Already initialized
            return;

        } else {
            AEM = new ArchitectureEvaluationManager(params, evaluator);
            this.paramsMap.put(key, params);
            this.architectureEvaluationManagerMap.put(key, AEM);

            // Initialization
            AEM.init(1);
        }
    }

    private BaseParams getProblemParameters(String problem) {
        String key;

        if (problem.equalsIgnoreCase("SMAP")) {
            key = "SMAP";
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

        if (problem.equalsIgnoreCase("SMAP") || problem.equalsIgnoreCase("ClimateCentric")) {
            // Generate a new architecture
            architecture = new rbsa.eoss.problems.Assigning.Architecture(bitString, numSatellites, (rbsa.eoss.problems.Assigning.AssigningParams) params);

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
            architecture = new rbsa.eoss.problems.PartitioningAndAssigning.Architecture(instrPartitioning, orbitAssignation, numSatellites, (rbsa.eoss.problems.PartitioningAndAssigning.Decadal2017AerosolsParams) params);

        } else {
            throw new IllegalArgumentException("Unrecorgnizable problem type: " + problem);
        }

        return architecture;
    }

    @Override
    public BinaryInputArchitecture evalBinaryInputArch(String problem, List<Boolean> boolList) {

        // Initialize Jess
        initJess(problem);

        // Input a new architecture design
        String bitString = "";
        for (Boolean b : boolList) {
            bitString += b ? "1" : "0";
        }

        BaseParams params = this.paramsMap.get(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        AbstractArchitecture architecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);

        // Evaluate the architecture
        Result result = AEM.evaluateArchitecture(architecture, "Slow");

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

        // Initialize Jess
        initJess(problem);

        // Input a new architecture design
        int[] intArray = new int[intList.size()];
        for (int i = 0; i < intList.size(); i++) {
            intArray[i] = intList.get(i);
        }

        BaseParams params = this.paramsMap.get(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);
        AbstractArchitecture architecture = this.getArchitectureDiscreteInput(problem, intArray, 1, params);

        // Evaluate the architecture
        Result result = AEM.evaluateArchitecture(architecture, "Slow");

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
    public List<BinaryInputArchitecture> runLocalSearchBinaryInput(String problem, List<Boolean> boolList) {

        BaseParams params = this.paramsMap.get(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);

        String bitString = "";
        for (Boolean b : boolList) {
            bitString += b ? "1" : "0";
        }

        ArrayList<String> samples = randomLocalChangeBinaryInput(bitString, 4, params);

        List<BinaryInputArchitecture> out = new ArrayList<>();
        for (String sample : samples) {

            AbstractArchitecture architecture = this.getArchitectureBinaryInput(problem, sample, 1, params);

            // Evaluate the architecture
            Result result = AEM.evaluateArchitecture(architecture, "Slow");

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
    public List<DiscreteInputArchitecture> runLocalSearchDiscreteInput(String problem, List<Integer> boolList) {
        throw new UnsupportedOperationException("Local search for discrete input is not supported yet.");
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
    public List<String> getCritiqueBinaryInputArch(String problem, List<Boolean> boolList) {

        String bitString = "";
        for (Boolean b : boolList) {
            bitString += b ? "1" : "0";
        }

        AbstractArchitecture architecture;

        System.out.println(bitString);

        if (problem.equalsIgnoreCase("SMAP") || problem.equalsIgnoreCase("ClimateCentric")) {

            rbsa.eoss.problems.Assigning.AssigningParams params = (rbsa.eoss.problems.Assigning.AssigningParams) this.paramsMap.get(problem);
            ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);

            // Generate a new architecture
            architecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);

            // Initialize Critique Generator
            rbsa.eoss.problems.Assigning.CritiqueGenerator critiquer = new rbsa.eoss.problems.Assigning.CritiqueGenerator(params, AEM.getResourcePool(), architecture);

            return critiquer.getCritique();

        } else {
            throw new IllegalArgumentException("Unrecorgnizable problem type: " + problem);
        }
    }

    @Override
    public List<String> getCritiqueDiscreteInputArch(String problem, List<Integer> boolList) {
        throw new UnsupportedOperationException("getCritiqueDiscreteArch() not supported yet.");
    }

    @Override
    public ArrayList<String> getOrbitList(String problem) {
        BaseParams params = this.getProblemParameters(problem);
        ArrayList<String> orbitList = new ArrayList<>();
        for (String o : params.getOrbitList()) {
            orbitList.add(o);
        }
        return orbitList;
    }

    @Override
    public ArrayList<String> getInstrumentList(String problem) {
        BaseParams params = this.getProblemParameters(problem);
        ArrayList<String> instrumentList = new ArrayList<>();
        for (String i : params.getInstrumentList()) {
            instrumentList.add(i);
        }
        return instrumentList;
    }

    @Override
    public ArrayList<String> getObjectiveList(String problem) {
        BaseParams params = this.getProblemParameters(problem);
        ArrayList<String> objectiveList = new ArrayList<>();
        params.objectiveDescriptions.forEach((k, v) -> {
            objectiveList.add(k);
        });
        return objectiveList;
    }

    @Override
    public ArrayList<String> getInstrumentsForObjective(String problem, String objective) {
        BaseParams params = this.getProblemParameters(problem);
        return new ArrayList<>(params.objectivesToInstruments.get(objective));
    }

    @Override
    public ArrayList<String> getInstrumentsForPanel(String problem, String panel) {
        BaseParams params = this.getProblemParameters(problem);
        return new ArrayList<>(params.panelsToInstruments.get(panel));
    }

    @Override
    public List<ObjectiveSatisfaction> getArchitectureScoreExplanation(String problem, List<Boolean> arch) {
        String bitString = "";
        for (Boolean b : arch) {
            bitString += b ? "1" : "0";
        }

        BaseParams params = this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);

        // Generate a new architecture
        AbstractArchitecture architecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);

        // Evaluate the architecture
        Result result = null;
        // Save the explanations for each stakeholder score
        List<ObjectiveSatisfaction> explanations = new ArrayList<>();

        result = AEM.evaluateArchitecture(architecture, "Slow");
        for (int i = 0; i < params.panelNames.size(); ++i) {
            explanations.add(new ObjectiveSatisfaction(params.panelNames.get(i),
                    result.getPanelScores().get(i), params.panelWeights.get(i)));
        }

        return explanations;
    }


    @Override
    public List<ObjectiveSatisfaction> getPanelScoreExplanation(String problem, List<Boolean> arch, String panel) {
        String bitString = "";
        for (Boolean b : arch) {
            bitString += b ? "1" : "0";
        }

        BaseParams params = this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);

        // Generate a new architecture
        AbstractArchitecture architecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);

        // Evaluate the architecture
        Result result = null;
        // Save the explanations for each stakeholder score
        List<ObjectiveSatisfaction> explanations = new ArrayList<>();

        result = AEM.evaluateArchitecture(architecture, "Slow");
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
    public List<ObjectiveSatisfaction> getObjectiveScoreExplanation(String problem, List<Boolean> arch, String objective) {
        String bitString = "";
        for (Boolean b : arch) {
            bitString += b ? "1" : "0";
        }

        BaseParams params = this.getProblemParameters(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);

        // Generate a new architecture
        AbstractArchitecture architecture = this.getArchitectureBinaryInput(problem, bitString, 1, params);

        // Evaluate the architecture
        Result result = null;
        // Save the explanations for each stakeholder score
        List<ObjectiveSatisfaction> explanations = new ArrayList<>();

        result = AEM.evaluateArchitecture(architecture, "Slow");
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

    @Override
    public SubobjectiveDetails getSubscoreDetailsBinaryInput(String problem, BinaryInputArchitecture architecture, String subobj) {
        throw new UnsupportedOperationException();
    }

    @Override
    public SubobjectiveDetails getSubscoreDetailsDiscreteInput(String problem, DiscreteInputArchitecture architecture, String subobj) {
        throw new UnsupportedOperationException();
    }

    @Override
    public List<SubscoreInformation> getArchScienceInformationBinaryInput(String problem, BinaryInputArchitecture architecture) {
        throw new UnsupportedOperationException();
    }

    @Override
    public List<SubscoreInformation> getArchScienceInformationDiscreteInput(String problem, DiscreteInputArchitecture architecture) {
        throw new UnsupportedOperationException();
    }

    @Override
    public List<MissionCostInformation> getArchCostInformationBinaryInput(String problem, BinaryInputArchitecture architecture) {
        throw new UnsupportedOperationException();
    }

    @Override
    public List<MissionCostInformation> getArchCostInformationDiscreteInput(String problem, DiscreteInputArchitecture architecture) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void startGABinaryInput(String problem, List<BinaryInputArchitecture> dataset, String username) {
        //PATH
        String path = ".";

        ExecutorService pool = Executors.newFixedThreadPool(8);
        CompletionService<Algorithm> ecs = new ExecutorCompletionService<>(pool);

        //parameters and operators for search
        TypedProperties properties = new TypedProperties();
        //search paramaters set here
        int popSize = 10;
        int maxEvals = 50;
        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);
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

        //setup for saving results
        properties.setBoolean("saveQuality", true);
        properties.setBoolean("saveCredits", true);
        properties.setBoolean("saveSelection", true);

        //initialize problem
        BaseParams params = this.paramsMap.get(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);

        AEM.reset();
        AEM.init(8);
        Problem assignmentProblem = new AssigningProblem(new int[]{1}, problem, AEM, params);

        // Create a solution for each input arch in the dataset
        List<Solution> initial = new ArrayList<>(dataset.size());
        for (int i = 0; i < dataset.size(); ++i) {
            AssigningArchitecture new_arch = new AssigningArchitecture(new int[]{1},
                    params.getNumInstr(), params.getNumOrbits(), 2);
            for (int j = 1; j < new_arch.getNumberOfVariables(); ++j) {
                BinaryVariable var = new BinaryVariable(1);
                var.set(0, dataset.get(i).inputs.get(j));
                new_arch.setVariable(j, var);
            }
            new_arch.setObjective(0, dataset.get(i).outputs.get(0));
            new_arch.setObjective(1, dataset.get(i).outputs.get(1));
            initial.set(i, new_arch);
        }
        initialization = new InjectedInitialization(assignmentProblem, popSize, initial);

        //initialize population structure for algorithm
        Population population = new Population();
        EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
        ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
        TournamentSelection selection = new TournamentSelection(2, comp);

        singlecross = new OnePointCrossover(crossoverProbability);
        bitFlip = new BitFlip(mutationProbability);
        intergerMutation = new IntegerUM(mutationProbability);
        CompoundVariation var = new CompoundVariation(singlecross, bitFlip, intergerMutation);

        // REDIS
        RedisClient redisClient = RedisClient.create("redis://localhost:6379/0");

        Algorithm eMOEA = new EpsilonMOEA(assignmentProblem, population, archive, selection, var, initialization);
        ecs.submit(new InteractiveSearch(eMOEA, properties, username, redisClient));

        try {
            Algorithm alg = ecs.take().get();
        } catch (InterruptedException | ExecutionException ex) {
            ex.printStackTrace();
        }

        redisClient.shutdown();
        AEM.clear();
        pool.shutdown();
        System.out.println("DONE");
    }

    @Override
    public void startGADiscreteInput(String problem, List<DiscreteInputArchitecture> dataset, String username) {
        //PATH
        String path = ".";

        ExecutorService pool = Executors.newFixedThreadPool(8);
        CompletionService<Algorithm> ecs = new ExecutorCompletionService<>(pool);

        //parameters and operators for search
        TypedProperties properties = new TypedProperties();
        //search paramaters set here
        int popSize = 10;
        int maxEvals = 50;
        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);

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

        //setup for saving results
        properties.setBoolean("saveQuality", true);
        properties.setBoolean("saveCredits", true);
        properties.setBoolean("saveSelection", true);

        //initialize problem
        BaseParams params = this.paramsMap.get(problem);
        ArchitectureEvaluationManager AEM = this.architectureEvaluationManagerMap.get(problem);

        AEM.reset();
        AEM.init(8);
        Problem partitioningAndAssigningProblem = new search.problems.PartitioningAndAssigning.PartitioningAndAssigningProblem(problem, AEM, params);

        // Create a solution for each input arch in the dataset
        List<Solution> initial = new ArrayList<>(dataset.size());
        for (int i = 0; i < dataset.size(); ++i) {
            PartitioningAndAssigningArchitecture new_arch = new PartitioningAndAssigningArchitecture(
                    params.getNumInstr(), params.getNumOrbits(), 2);

            int numPartitioningVariables = new_arch.getDecision("instrumentPartitioning").getNumberOfVariables();
            int numAssignmentVariables = new_arch.getDecision("orbitAssignment").getNumberOfVariables();

            for (int j = 0; j < numPartitioningVariables; ++j) {
                IntegerVariable var = new IntegerVariable(dataset.get(i).inputs.get(j), 0, params.getNumInstr());
                new_arch.setVariable(j, var);
            }

            for (int j = numPartitioningVariables; j < numPartitioningVariables + numAssignmentVariables; ++j) {
                IntegerVariable var = new IntegerVariable(dataset.get(i).inputs.get(j), 0, params.getNumOrbits());
                new_arch.setVariable(j, var);
            }

            new_arch.setObjective(0, dataset.get(i).outputs.get(0));
            new_arch.setObjective(1, dataset.get(i).outputs.get(1));
            initial.set(i, new_arch);
        }

        initialization = new InjectedInitialization(partitioningAndAssigningProblem, popSize, initial);

        //initialize population structure for algorithm
        Population population = new Population();
        EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
        ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
        TournamentSelection selection = new TournamentSelection(2, comp);

        singlecross = new OnePointCrossover(crossoverProbability);
        bitFlip = new BitFlip(mutationProbability);
        intergerMutation = new IntegerUM(mutationProbability);
        CompoundVariation var = new CompoundVariation(singlecross, bitFlip, intergerMutation);

        // REDIS
        RedisClient redisClient = RedisClient.create("redis://localhost:6379/0");

        Algorithm eMOEA = new EpsilonMOEA(partitioningAndAssigningProblem, population, archive, selection, var, initialization);
        ecs.submit(new InteractiveSearch(eMOEA, properties, username, redisClient));

        try {
            Algorithm alg = ecs.take().get();
        } catch (InterruptedException | ExecutionException ex) {
            ex.printStackTrace();
        }

        redisClient.shutdown();
        AEM.clear();
        pool.shutdown();
        System.out.println("DONE");
    }
}
