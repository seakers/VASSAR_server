package test;

import org.moeaframework.core.Initialization;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.operator.CompoundVariation;
import org.moeaframework.util.TypedProperties;
import rbsa.eoss.evaluation.AbstractArchitectureEvaluator;
import rbsa.eoss.evaluation.ArchitectureEvaluationManager;
import rbsa.eoss.local.BaseParams;
import seak.architecture.Architecture;
import seak.architecture.util.IntegerVariable;
import search.problems.PartitioningAndAssigning.PartitioningAndAssigningArchitecture;
import search.problems.PartitioningAndAssigning.PartitioningAndAssigningInitialization;

import java.io.File;
import java.util.*;


public class RepairTest {

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

        Problem partitioningAndAssigningProblem = new search.problems.PartitioningAndAssigning.PartitioningAndAssigningProblem(problem, AEM, params);

        // Set params
        int popSize = 300;
        double mutationProp = 1. / 7.;

        Initialization initialization = new PartitioningAndAssigningInitialization(partitioningAndAssigningProblem, popSize, params);

        // Generate initial solutions
        Solution[] testSolutions = initialization.initialize();

        Variation mutation = new search.problems.PartitioningAndAssigning.operators.PartitioningAndAssigningMutation(mutationProp, params);
        Variation operator = mutation;

        int testNum = 100;
        Random random = new Random();

        for(int i = 0; i < testNum; i++){

            Solution[] parents;
            Solution[] children;

            int a = random.nextInt(testSolutions.length);
            parents = new Solution[1];
            parents[0] = testSolutions[a];

            children = operator.evolve(parents);
            int[] inputs = getIntVariables((search.problems.PartitioningAndAssigning.PartitioningAndAssigningArchitecture) children[0]);

            if(!isFeasibleAssignment(inputs, params)){
                String temp = ppArchInput(inputs, params);
                System.out.println("Before repair:  " + temp);

                int[] repairedInputs = repair(inputs, params);
                String temp2 = ppArchInput(repairedInputs, params);
                System.out.println("After  repair:  " + temp2);
                System.out.println("----");
            }
        }

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

    public static boolean isFeasibleAssignment(int[] inputs, BaseParams params) {

        int[] partitioning = new int[params.getNumInstr()];
        int[] assignment = new int[params.getNumInstr()];
        for(int i = 0; i < inputs.length; i++){
            if(i < params.getNumInstr()){
                partitioning[i] = inputs[i];
            }else{
                assignment[i - params.getNumInstr()] = inputs[i];
            }
        }

        Set<Integer> temp = new HashSet<>();
        for(int p: partitioning){
            if(p == -1){
                // Instrument not assigned to any set
                return false;
            }
            temp.add(p);
        }

        int numOrbitUsed = 0;
        for(int o: assignment){
            if(o!=-1){
                numOrbitUsed++;
            }
        }
        if(numOrbitUsed != temp.size()){
            // Number of orbits assigned does not match the number of satellites
            return false;
        }

        for(int p:temp){
            if(assignment[p] == -1){
                // Set not assigned to any orbit
                return false;
            }
        }
        return true;
    }

    public static int[] repair(int[] inputs, BaseParams params){

        int[] partitioning = Arrays.copyOfRange(inputs, 0, params.getNumInstr());
        int[] assigning = Arrays.copyOfRange(inputs, params.getNumInstr(), 2 * params.getNumInstr()+1);

        Architecture newArch = new PartitioningAndAssigningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2);
        int[] newPartitioning = new int[partitioning.length];
        int[] newAssigning = new int[assigning.length];

        HashMap<Integer, Integer> satID2satIndex = new HashMap<>();
        HashMap<Integer, Integer> satIndex2Orbit = new HashMap<>();
        int satIndex = 0;
        for(int i = 0; i < partitioning.length; i++){
            int satID = partitioning[i];
            if(!satID2satIndex.containsKey(satID)){
                satID2satIndex.put(satID, satIndex);
                satIndex2Orbit.put(satIndex, assigning[satID]);
                satIndex++;
            }

            newPartitioning[i] = satID2satIndex.get(satID);
        }

        for(int i = 0; i < assigning.length;i ++){

            if(satIndex2Orbit.containsKey(i)){
                int orb = satIndex2Orbit.get(i);
                if(orb == -1){
                    Collection<Integer> orbitsUsed = satIndex2Orbit.values();
                    ArrayList<Integer> orbitOptions;

                    if(orbitsUsed.size() == params.getNumOrbits()){
                        orbitOptions = new ArrayList<>(orbitsUsed);

                    }else{
                        orbitOptions = new ArrayList<>();
                        for(int j = 0; j < params.getNumOrbits(); j++){
                            if(!orbitsUsed.contains(j)){
                                orbitOptions.add(j);
                            }
                        }
                    }

                    Collections.shuffle(orbitOptions);
                    orb = orbitOptions.get(0);
                    satIndex2Orbit.put(i, orb);
                }
                newAssigning[i] = orb;

            }else{
                newAssigning[i] = -1;
            }
        }

        int[] newIntVars = new int[partitioning.length + assigning.length];
        for(int i = 0; i < partitioning.length;i ++){
            newIntVars[i] = newPartitioning[i];
        }
        for(int i = 0; i < assigning.length;i ++){
            newIntVars[i + newPartitioning.length] = newAssigning[i];
        }
        return newIntVars;
    }
}
