package test;

import io.lettuce.core.RedisClient;
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
import search.BinaryInputInteractiveSearch;
import search.problems.PartitioningAndAssigning.PartitioningAndAssigningInitialization;

import java.io.File;


public class Test {

    public static void main(String[] args){

        String problem = "Decadal2017Aerosols";
        String search_clps = "";
        String root = System.getProperty("user.dir");

        //PATH
        String path = root +
                File.separator + "problems" +
                File.separator + "SMAP";

        BaseParams params = new rbsa.eoss.problems.PartitioningAndAssigning.Decadal2017AerosolsParams(path, "CRISP-ATTRIBUTES", "test", "normal", search_clps);
        AbstractArchitectureEvaluator evaluator = new rbsa.eoss.problems.PartitioningAndAssigning.ArchitectureEvaluator(params);
        ArchitectureEvaluationManager AEM = new ArchitectureEvaluationManager(params, evaluator);

        //parameters and operators for search
        TypedProperties properties = new TypedProperties();
        //search paramaters set here
        int popSize = 20;
        int maxEvals = 50;

        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);

        double crossoverProbability = 1.0;
        properties.setDouble("crossoverProbability", crossoverProbability);
        double mutationProbability = 1. / 60.;
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
        Problem partitioningAndAssigningProblem = new search.problems.PartitioningAndAssigning.PartitioningAndAssigningProblem(problem, AEM, params);


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
        Variation singlecross = new search.problems.PartitioningAndAssigning.operators.PartitioningAndAssigningCrossover(crossoverProbability, params);
        Variation intergerMutation = new search.problems.PartitioningAndAssigning.operators.PartitioningAndAssigningMutation(mutationProbability, params);
        CompoundVariation var = new CompoundVariation(singlecross, intergerMutation);

        // REDIS
        RedisClient redisClient = RedisClient.create("redis://localhost:6379/0");

        Algorithm eMOEA = new EpsilonMOEA(partitioningAndAssigningProblem, population, archive, selection, var, initialization);
        new BinaryInputInteractiveSearch(eMOEA, properties, "hbang", redisClient).call();

        redisClient.shutdown();
        AEM.clear();
        System.out.println("DONE");
    }
}
