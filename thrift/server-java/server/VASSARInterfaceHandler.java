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


import java.util.*;
import java.util.concurrent.*;

import io.lettuce.core.RedisClient;
import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.util.TypedProperties;
import rbsa.eoss.*;
import rbsa.eoss.local.Params;
import javaInterface.BinaryInputArchitecture;
import javaInterface.VASSARInterface;
import javaInterface.ObjectiveSatisfaction;
import seak.architecture.operators.IntegerUM;

public class VASSARInterfaceHandler implements VASSARInterface.Iface {

    private Params params;
    private ArchitectureEvaluator AE = null;

    public VASSARInterfaceHandler() {
        initJess();
    }

    public void ping() {
      System.out.println("ping()");
    }
  
    private void initJess() {
        // Set a path to the project folder
        String path = System.getProperty("user.dir");
        
        // Initialization
        String search_clps = "";
        params = Params.initInstance(path, "FUZZY-ATTRIBUTES", "test","normal", search_clps);//FUZZY or CRISP
        AE = ArchitectureEvaluator.getInstance();
        AE.init(1);
    }

    @Override
    public BinaryInputArchitecture eval(List<Boolean> boolList) {
        // Input a new architecture design
        // There must be 5 orbits. Instrument name is represented by a capital letter, taken from {A,B,C,D,E,F,G,H,I,J,K,L}
        
        String bitString = "";
        for (Boolean b: boolList) {
            bitString += b ? "1" : "0";
        }

        // Generate a new architecture
        Architecture architecture = new Architecture(bitString, 1);

        // Evaluate the architecture
        Result result = AE.evaluateArchitecture(architecture,"Slow");
        
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
    public List<BinaryInputArchitecture> runLocalSearch(List<Boolean> boolList) {
        String bitString = "";
        for (Boolean b: boolList) {
            bitString += b ? "1" : "0";
        }

        ArrayList<String> samples = randomLocalChange(bitString, 4);

        List<BinaryInputArchitecture> out = new ArrayList<>();

        for (String sample: samples) {
            // Generate a new architecture
            Architecture architecture = new Architecture(sample, 1);

            // Evaluate the architecture
            Result result = AE.evaluateArchitecture(architecture,"Slow");

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

    private ArrayList<String> randomLocalChange(String bitString, int n) {
        Random rand = new Random();
        int numVars = params.orbitList.length * params.instrumentList.length;

        ArrayList<String> out = new ArrayList<>();

        for (int i = 0; i < n; i++) {
            int k = rand.nextInt(numVars);

            StringBuilder tempBitString = new StringBuilder(bitString);
            if (bitString.charAt(k) == '1') {
                tempBitString.setCharAt(k, '0');
            }
            else {
                tempBitString.setCharAt(k, '1');
            }
            out.add(tempBitString.toString());
        }
        return out;
    }

    private List<Boolean> bitString2BoolArray(String bitString){
        List<Boolean> out = new ArrayList<>();
        for (int i = 0; i < bitString.length(); i++) {
            out.add(bitString.charAt(i) == '1');
        }
        return out;
    }

    @Override
    public List<String> getCritique(List<Boolean> boolList) {
        String bitString = "";
        for(Boolean b: boolList){
            bitString += b ? "1" : "0";
        }
        
        System.out.println(bitString);

        // Generate a new architecture
        Architecture architecture = new Architecture(bitString, 1);

        // Initialize Critique Generator
        CritiqueGenerator critiquer = new CritiqueGenerator(architecture);

        return critiquer.getCritique();
    }

    @Override
    public ArrayList<String> getOrbitList() {
        ArrayList<String> orbitList = new ArrayList<>();
        for(String o: params.orbitList){
            orbitList.add(o);
        }
        return orbitList;
    }

    @Override
    public ArrayList<String> getInstrumentList() {
        ArrayList<String> instrumentList = new ArrayList<>();
        for (String i: params.instrumentList) {
            instrumentList.add(i);
        }
        return instrumentList;
    }

    @Override
    public ArrayList<String> getObjectiveList() {
        ArrayList<String> objectiveList = new ArrayList<>();
        params.objectiveDescriptions.forEach((k, v) -> {
            objectiveList.add(k);
        });
        return objectiveList;
    }

    @Override
    public ArrayList<String> getInstrumentsForObjective(String objective) {
        return new ArrayList<>(params.objectivesToInstruments.get(objective));
    }

    @Override
    public ArrayList<String> getInstrumentsForPanel(String panel) {
        return new ArrayList<>(params.panelsToInstruments.get(panel));
    }

    @Override
    public List<ObjectiveSatisfaction> getArchitectureScoreExplanation(List<Boolean> arch) {
        String bitString = "";
        for (Boolean b: arch) {
            bitString += b ? "1" : "0";
        }

        // Generate a new architecture
        Architecture architecture = new Architecture(bitString, 1);
        architecture.setEvalMode("DEBUG");

        // Evaluate the architecture
        Result result = null;
        // Save the explanations for each stakeholder score
        List<ObjectiveSatisfaction> explanations = new ArrayList<>();

        result = AE.evaluateArchitecture(architecture, "Slow");
        for (int i = 0; i < params.panelNames.size(); ++i) {
            explanations.add(new ObjectiveSatisfaction(params.panelNames.get(i),
                    result.getPanelScores().get(i), params.panelWeights.get(i)));
        }

        return explanations;
    }


    @Override
    public List<ObjectiveSatisfaction> getPanelScoreExplanation(List<Boolean> arch, String panel) {
        String bitString = "";
        for (Boolean b: arch) {
            bitString += b ? "1" : "0";
        }

        // Generate a new architecture
        Architecture architecture = new Architecture(bitString, 1);
        architecture.setEvalMode("DEBUG");

        // Evaluate the architecture
        Result result = null;
        // Save the explanations for each stakeholder score
        List<ObjectiveSatisfaction> explanations = new ArrayList<>();

        result = AE.evaluateArchitecture(architecture, "Slow");
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
    public List<ObjectiveSatisfaction> getObjectiveScoreExplanation(List<Boolean> arch, String objective) {
        String bitString = "";
        for (Boolean b: arch) {
            bitString += b ? "1" : "0";
        }

        // Generate a new architecture
        Architecture architecture = new Architecture(bitString, 1);
        architecture.setEvalMode("DEBUG");

        // Evaluate the architecture
        Result result = null;
        // Save the explanations for each stakeholder score
        List<ObjectiveSatisfaction> explanations = new ArrayList<>();

        result = AE.evaluateArchitecture(architecture, "Slow");
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
    public void startGA(List<BinaryInputArchitecture> dataset, String username) {
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
        Params.initInstance(path, "CRISP-ATTRIBUTES", "test","normal","");
        ArchitectureEvaluator.getInstance().init(8);
        Problem problem = new InstrumentAssignment(new int[]{1});

        // Create a solution for each input arch in the dataset
        List<Solution> initial = new ArrayList<>(dataset.size());
        for (int i = 0; i < dataset.size(); ++i) {
            InstrumentAssignmentArchitecture new_arch = new InstrumentAssignmentArchitecture(new int[]{1},
                    Params.getInstance().numInstr, Params.getInstance().numOrbits, 2);
            for (int j = 1; j < new_arch.getNumberOfVariables(); ++j) {
                BinaryVariable var = new BinaryVariable(1);
                var.set(0,dataset.get(i).inputs.get(j));
                new_arch.setVariable(j, var);
            }
            new_arch.setObjective(0, dataset.get(i).outputs.get(0));
            new_arch.setObjective(1, dataset.get(i).outputs.get(1));
            initial.set(i, new_arch);
        }
        initialization = new InjectedInitialization(problem, popSize, initial);

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

        Algorithm eMOEA = new EpsilonMOEA(problem, population, archive, selection, var, initialization);
        ecs.submit(new InteractiveSearch(eMOEA, properties, username, redisClient));

        try {
            Algorithm alg = ecs.take().get();
        } catch (InterruptedException | ExecutionException ex) {
            ex.printStackTrace();
        }

        redisClient.shutdown();
        ArchitectureEvaluator.getInstance().clear();
        pool.shutdown();
        System.out.println("DONE");
    }
}

