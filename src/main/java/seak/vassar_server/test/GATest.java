package seak.vassar_server.test;

import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.CompoundVariation;
import org.moeaframework.core.operator.TournamentSelection;
import org.moeaframework.util.TypedProperties;
import rbsa.eoss.evaluation.AbstractArchitectureEvaluator;
import rbsa.eoss.evaluation.ArchitectureEvaluationManager;
import rbsa.eoss.local.BaseParams;
import seak.architecture.util.IntegerVariable;
import seak.vassar_server.search.problems.PartitioningAndAssigning.PartitioningAndAssigningArchitecture;
import seak.vassar_server.search.problems.PartitioningAndAssigning.PartitioningAndAssigningInitialization;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.StringJoiner;


public class GATest {

    public static void main(String[] args){

        String problem = "Decadal2017Aerosols";
        String search_clps = "";
        String root = System.getProperty("user.dir");

        //PATH
        String path = root +
                File.separator + "problems" +
                File.separator + "SMAP";

        BaseParams params = new rbsa.eoss.problems.PartitioningAndAssigning.Decadal2017AerosolsParams(path, "CRISP-ATTRIBUTES", "seak/vassar_server/test", "normal", search_clps);
        AbstractArchitectureEvaluator evaluator = new rbsa.eoss.problems.PartitioningAndAssigning.ArchitectureEvaluator(params);
        ArchitectureEvaluationManager AEM = new ArchitectureEvaluationManager(params, evaluator);

        //parameters and operators for seak.vassar_server.search
        TypedProperties properties = new TypedProperties();
        //seak.vassar_server.search paramaters set here
        int popSize = 100;
        int maxEvals = 1400;

        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);

        double crossoverProbability = 1.0;
        properties.setDouble("crossoverProbability", crossoverProbability);
        double mutationProbability = 1. / 10.;
        properties.setDouble("mutationProbability", mutationProbability);

        //setup for epsilon MOEA
        double[] epsilonDouble = new double[]{0.001, 1};

        //setup for saving results
        properties.setBoolean("saveQuality", true);
        properties.setBoolean("saveCredits", true);
        properties.setBoolean("saveSelection", true);

        //initialize problem

        AEM.reset();
        AEM.init(1);
        Problem partitioningAndAssigningProblem = new seak.vassar_server.search.problems.PartitioningAndAssigning.PartitioningAndAssigningProblem(problem, AEM, params);


//        // Create a solution for each input arch in the dataset
//        List<Solution> initial = new ArrayList<>(dataset.size());
//        for (int i = 0; i < dataset.size(); ++i) {
//            PartitioningAndAssigningArchitecture new_arch = new PartitioningAndAssigningArchitecture(
//                    params.getNumInstr(), params.getNumOrbits(), 2);
//
//            int numPartitioningVariables = new_arch.getDecision("instrumentPartitioning").getNumberOfVariables();
//            int numAssignmentVariables = new_arch.getDecision("orbitAssignment").getNumberOfVariables();
//
//            for (int j = 0; j < numPartitioningVariables; ++j) {
//                IntegerVariable var = new IntegerVariable(dataset.get(i).inputs.get(j), 0, params.getNumInstr());
//                new_arch.setVariable(j, var);
//            }
//
//            for (int j = numPartitioningVariables; j < numPartitioningVariables + numAssignmentVariables; ++j) {
//                IntegerVariable var = new IntegerVariable(dataset.get(i).inputs.get(j), -1, params.getNumOrbits());
//                new_arch.setVariable(j, var);
//            }
//
//            new_arch.setObjective(0, dataset.get(i).outputs.get(0));
//            new_arch.setObjective(1, dataset.get(i).outputs.get(1));
//            initial.add(new_arch);
//        }

        Initialization initialization = new PartitioningAndAssigningInitialization(partitioningAndAssigningProblem, popSize, params);

        //initialize population structure for algorithm
        Population population = new Population();
        EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
        ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
        TournamentSelection selection = new TournamentSelection(2, comp);

        // Define operators
        Variation singlecross = new seak.vassar_server.search.problems.PartitioningAndAssigning.operators.PartitioningAndAssigningCrossover(crossoverProbability, params);
        Variation intergerMutation = new seak.vassar_server.search.problems.PartitioningAndAssigning.operators.PartitioningAndAssigningMutation(mutationProbability, params);
        CompoundVariation var = new CompoundVariation(singlecross, intergerMutation);

        Algorithm eMOEA = new EpsilonMOEA(partitioningAndAssigningProblem, population, archive, selection, var, initialization);

        Algorithm alg = eMOEA;
        int populationSize = popSize;
        int maxEvaluations = maxEvals;

        // run the executor using the listener to collect results
        System.out.println("Starting " + alg.getClass().getSimpleName() + " on " + alg.getProblem().getName() + " with pop size: " + populationSize);
        alg.step();
        long startTime = System.currentTimeMillis();

        while (!alg.isTerminated() && (alg.getNumberOfEvaluations() < maxEvaluations)) {
            alg.step();
        }

        Population pop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();

        int numInputs = 10;
        ArrayList<ArrayList<Integer>> inputs = new ArrayList<>();
        ArrayList<ArrayList<Double>> outputs = new ArrayList<>();
        for(int i = 0; i < pop.size(); i++){
            PartitioningAndAssigningArchitecture sol = (PartitioningAndAssigningArchitecture)pop.get(i);
            numInputs = sol.getNumberOfVariables();
            ArrayList<Integer> input = new ArrayList<>();
            ArrayList<Double> output = new ArrayList<>();
            for(int j = 0; j < sol.getNumberOfVariables(); j++){
                input.add(((IntegerVariable) sol.getVariable(j)).getValue());
            }
            for(int j = 0; j < sol.getNumberOfObjectives(); j++){
                output.add(sol.getObjective(j));
            }
            inputs.add(input);
            outputs.add(output);
        }

        String filePath = path + "/results/GATest.csv";

        int count = 0;
        try (BufferedWriter outputWriter = new BufferedWriter(new FileWriter(filePath))) {

            StringJoiner header = new StringJoiner(",");
            for(int i = 0; i < numInputs; i++){
                header.add("input" + i);
            }

            header.add("Science");
            header.add("Cost");
            outputWriter.write(header.toString() + "\n");

            for (int i = 0; i < inputs.size(); i++) {
                if (outputs.get(i).get(0) == 0.0) {
                    count++;

                }else {
                    StringJoiner sj = new StringJoiner(",");
                    StringJoiner inputString = new StringJoiner(",");
                    for(int j:inputs.get(i)){
                        String temp = "" + j;
                        inputString.add(temp);
                    }
                    sj.add(inputString.toString());
                    sj.add(Double.toString(outputs.get(i).get(0)));
                    sj.add(Double.toString(outputs.get(i).get(1)));
                    outputWriter.write(sj.toString() + "\n");
                }
            }
        }catch(Exception e){
            e.printStackTrace();
        }


        alg.terminate();
        long finishTime = System.currentTimeMillis();
        System.out.println("Done with optimization. Execution time: " + ((finishTime - startTime) / 1000) + "s");

        AEM.clear();
        System.out.println("DONE");
    }
}
