package seakers.engineerserver.test;

import org.moeaframework.core.*;
import org.moeaframework.core.operator.CompoundVariation;
import org.moeaframework.util.TypedProperties;
import seakers.architecture.Architecture;
import seakers.architecture.util.IntegerVariable;
import seakers.engineerserver.search.problems.PartitioningAndAssigning.PartitioningAndAssigningInitialization;
import seakers.engineerserver.search.problems.PartitioningAndAssigning.PartitioningAndAssigningArchitecture;
import seakers.engineerserver.search.problems.PartitioningAndAssigning.PartitioningAndAssigningProblem;
import seakers.engineerserver.search.problems.PartitioningAndAssigning.operators.PartitioningAndAssigningCrossover;
import seakers.engineerserver.search.problems.PartitioningAndAssigning.operators.PartitioningAndAssigningMutation;
import seakers.vassar.evaluation.AbstractArchitectureEvaluator;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.BaseParams;
import seakers.vassar.problems.PartitioningAndAssigning.ArchitectureEvaluator;
import seakers.vassar.problems.PartitioningAndAssigning.Decadal2017AerosolsParams;

import java.io.File;
import java.util.Random;
import java.util.StringJoiner;


public class OperatorTest {

    public static void main(String[] args){

        String problem = "Decadal2017Aerosols";
        String search_clps = "";
        String root = System.getProperty("user.dir");

        //PATH
        String path = root +
                File.separator + "problems" +
                File.separator + "SMAP";

        BaseParams params = new Decadal2017AerosolsParams(path, "CRISP-ATTRIBUTES", "seakers/engineerserver/test", "normal");
        AbstractArchitectureEvaluator evaluator = new ArchitectureEvaluator();
        ArchitectureEvaluationManager AEM = new ArchitectureEvaluationManager(params, evaluator);

        //parameters and operators for search
        TypedProperties properties = new TypedProperties();
        //search paramaters set here

        Problem partitioningAndAssigningProblem = new PartitioningAndAssigningProblem(problem, AEM, params);

        // Set params
        int popSize = 300;
        double mutationProp = 1. / 7.;
        double crossoverProp = 1.0;

        Initialization initialization = new PartitioningAndAssigningInitialization(partitioningAndAssigningProblem, popSize, params);

        // Generate initial solutions
        Solution[] testSolutions = initialization.initialize();

        Variation mutation = new PartitioningAndAssigningMutation(mutationProp, params);
        Variation crossover = new PartitioningAndAssigningCrossover(crossoverProp, params);
        CompoundVariation var = new CompoundVariation(crossover, mutation);

        Variation operator = mutation;


//        for(int i = 0; i < testSolutions.length; i++){
//            Solution[] parents = new Solution[1];
//            parents[0] = testSolutions[i];
//
//            int[] inputs = getIntVariables((PartitioningAndAssigningArchitecture) parents[0]);
//            String temp = ppArchInput(inputs, params);
//            System.out.println("Parent: " + temp);
//
//            Solution[] children = mutation.evolve(parents);
//            int[] inputs2 = getIntVariables((PartitioningAndAssigningArchitecture) children[0]);
//            String temp2 = ppArchInput(inputs2, params);
//            System.out.println("Child:  " + temp2);
//
//            System.out.println("----");
//        }

        int testNum = 100;
        Random random = new Random();

        for(int i = 0; i < testNum; i++){

            Solution[] parents;
            Solution[] children;

            int a = random.nextInt(testSolutions.length);
            int b;
            if(operator.getArity() == 1){
                parents = new Solution[1];
                parents[0] = testSolutions[a];

                int[] inputs = getIntVariables((PartitioningAndAssigningArchitecture) parents[0]);
                String temp = ppArchInput(inputs, params);
                System.out.println("Parent: " + temp);

                children = operator.evolve(parents);
                inputs = getIntVariables((PartitioningAndAssigningArchitecture) children[0]);
                temp = ppArchInput(inputs, params);
                System.out.println("Child:  " + temp);

                System.out.println("----");

            }else{
                b = a;
                while(b == a){
                    b = random.nextInt(testSolutions.length);
                }
                parents = new Solution[2];

                parents[0] = testSolutions[a];
                parents[1] = testSolutions[b];

                int[] inputs = getIntVariables((PartitioningAndAssigningArchitecture) parents[0]);
                String temp = ppArchInput(inputs, params);
                System.out.println("Parent1: " + temp);

                inputs = getIntVariables((PartitioningAndAssigningArchitecture) parents[1]);
                temp = ppArchInput(inputs, params);
                System.out.println("Parent2: " + temp);

                children = operator.evolve(parents);
                inputs = getIntVariables((PartitioningAndAssigningArchitecture) children[0]);
                temp = ppArchInput(inputs, params);
                System.out.println("Child1:  " + temp);

                inputs = getIntVariables((PartitioningAndAssigningArchitecture) children[1]);
                temp = ppArchInput(inputs, params);
                System.out.println("Child2:  " + temp);

                System.out.println("----");
            }
        }


//
//        int numInputs = 10;
//        ArrayList<ArrayList<Integer>> inputs = new ArrayList<>();
//        ArrayList<ArrayList<Double>> outputs = new ArrayList<>();
//        for(int i = 0; i < pop.size(); i++){
//            PartitioningAndAssigningArchitecture sol = (PartitioningAndAssigningArchitecture)pop.get(i);
//            numInputs = sol.getNumberOfVariables();
//            ArrayList<Integer> input = new ArrayList<>();
//            ArrayList<Double> output = new ArrayList<>();
//            for(int j = 0; j < sol.getNumberOfVariables(); j++){
//                input.add(((IntegerVariable) sol.getVariable(j)).getValue());
//            }
//            for(int j = 0; j < sol.getNumberOfObjectives(); j++){
//                output.add(sol.getObjective(j));
//            }
//            inputs.add(input);
//            outputs.add(output);
//        }
//
//        String filePath = path + "/results/GATest.csv";
//
//        int count = 0;
//        try (BufferedWriter outputWriter = new BufferedWriter(new FileWriter(filePath))) {
//
//            StringJoiner header = new StringJoiner(",");
//            for(int i = 0; i < numInputs; i++){
//                header.add("input" + i);
//            }
//
//            header.add("Science");
//            header.add("Cost");
//            outputWriter.write(header.toString() + "\n");
//
//            for (int i = 0; i < inputs.size(); i++) {
//                if (outputs.get(i).get(0) == 0.0) {
//                    count++;
//
//                }else {
//                    StringJoiner sj = new StringJoiner(",");
//                    StringJoiner inputString = new StringJoiner(",");
//                    for(int j:inputs.get(i)){
//                        String temp = "" + j;
//                        inputString.add(temp);
//                    }
//                    sj.add(inputString.toString());
//                    sj.add(Double.toString(outputs.get(i).get(0)));
//                    sj.add(Double.toString(outputs.get(i).get(1)));
//                    outputWriter.write(sj.toString() + "\n");
//                }
//            }
//        }catch(Exception e){
//            e.printStackTrace();
//        }
//
//        alg.terminate();
//        long finishTime = System.currentTimeMillis();
//        System.out.println("Done with optimization. Execution time: " + ((finishTime - startTime) / 1000) + "s");

        //AEM.clear();
        System.out.println("DONE");
    }

    private static String ppArchInput(int[] inputs, BaseParams params){
        StringJoiner partitioning = new StringJoiner(",");
        StringJoiner assigning = new StringJoiner(",");
        for(int i = 0; i < params.getNumInstr(); i++){
            String temp = "" + inputs[i];
            partitioning.add(temp);
        }
        for(int i = 0; i < params.getNumInstr(); i++){
            String temp = "" + inputs[i + params.getNumInstr()];
            assigning.add(temp);
        }

        String out = partitioning.toString() + " | " + assigning.toString();
        return out;
    }

    private static int[] getIntVariables(Architecture arch){
        int[] out = new int[arch.getNumberOfVariables()];
        for(int i = 0; i < out.length; i++){
            out[i] = ((IntegerVariable) arch.getVariable(i)).getValue();
        }
        return out;
    }
}
