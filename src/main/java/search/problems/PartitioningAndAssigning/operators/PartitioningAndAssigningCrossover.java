package search.problems.PartitioningAndAssigning.operators;

import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.Variation;
import rbsa.eoss.local.BaseParams;
import seak.architecture.Architecture;
import seak.architecture.util.IntegerVariable;
import search.problems.PartitioningAndAssigning.PartitioningAndAssigningArchitecture;

import java.util.Random;
import java.util.HashMap;
import java.util.ArrayList;

public class PartitioningAndAssigningCrossover implements Variation{

    protected double probability;
    protected BaseParams params;
    protected String option = "antoni"; // antoni, prachi

    public PartitioningAndAssigningCrossover(double probability, BaseParams params, String option){
        this.probability = probability;
        this.params = params;
        this.option = option;
    }

    public PartitioningAndAssigningCrossover(double probability, BaseParams params){
        this.probability = probability;
        this.params = params;
    }

    @Override
    public Solution[] evolve(Solution[] parents){

        Solution[] out = new Solution[2];

        Architecture arch1 = (PartitioningAndAssigningArchitecture) parents[0];
        Architecture arch2 = (PartitioningAndAssigningArchitecture) parents[1];
        int[] partitioning1 = getIntVariables("instrumentPartitioning", arch1);
        int[] partitioning2 = getIntVariables("instrumentPartitioning", arch2);
        int[] assigning1 = getIntVariables("orbitAssignment", arch1);
        int[] assigning2 = getIntVariables("orbitAssignment", arch2);

        Architecture newArch1 = new PartitioningAndAssigningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2);
        Architecture newArch2 = new PartitioningAndAssigningArchitecture(params.getNumInstr(), params.getNumOrbits(), 2);
        int[] newPartitioning1 = new int[partitioning1.length];
        int[] newPartitioning2 = new int[partitioning2.length];
        int[] newAssigning1 = new int[assigning1.length];
        int[] newAssigning2 = new int[assigning2.length];

        Random random = new Random();
        int split = random.nextInt(partitioning1.length); // Single-point crossover

        if(option == "antoni"){
            int[] orbitAssigned1 = new int[partitioning1.length];
            for(int i = 0; i < partitioning1.length; i++){
                int sat = partitioning1[i];
                int orb = assigning1[sat];
                orbitAssigned1[i] = orb;
            }

            int[] orbitAssigned2 = new int[partitioning2.length];
            for(int i = 0; i < partitioning2.length; i++){
                int sat = partitioning2[i];
                int orb = assigning2[sat];
                orbitAssigned2[i] = orb;
            }

            ArrayList<int[]> orbitAssignedSwapped = swapSubarrays(orbitAssigned1, orbitAssigned2, split);

            int[] orbitAssignmentInfo1 = orbitAssignedSwapped.get(0);
            int[] orbitAssignmentInfo2 = orbitAssignedSwapped.get(1);

            ArrayList<int[]> temp1 = extractPartitioningAndAssigning(orbitAssignmentInfo1);
            ArrayList<int[]> temp2 = extractPartitioningAndAssigning(orbitAssignmentInfo2);

            newPartitioning1 = temp1.get(0);
            newAssigning1 = temp1.get(1);
            newPartitioning2 = temp2.get(0);
            newAssigning2 = temp2.get(0);

        }else if(option == "prachi"){
            throw new UnsupportedOperationException();
        }

        setIntVariables("instrumentPartitioning", newArch1, newPartitioning1);
        setIntVariables("orbitAssignment", newArch1, newAssigning1);
        setIntVariables("instrumentPartitioning", newArch2, newPartitioning2);
        setIntVariables("orbitAssignment", newArch2, newAssigning2);

        out[0] = newArch1;
        out[1] = newArch2;
        return out;
    }

    private ArrayList<int[]> extractPartitioningAndAssigning(int[] assignedOrbit){
        int[] partitioning = new int[assignedOrbit.length];
        int[] assigning = new int[assignedOrbit.length];

        for(int i = 0; i < assignedOrbit.length; i++){
            partitioning[i] = -1;
            assigning[i] = -1;
        }

        ArrayList<Integer> orbitUsed = new ArrayList<>();
        HashMap<Integer, Integer>  orbit2SatIndex = new HashMap<>();
        for(int i = 0; i < assignedOrbit.length; i++){
            int orb = assignedOrbit[i];
            if(!orbitUsed.contains(orb)){
                orbit2SatIndex.put(orb, orbitUsed.size());
                orbitUsed.add(orb);
            }

            partitioning[i] = orbit2SatIndex.get(orb);
        }

        for(int i = 0; i < orbitUsed.size(); i++){
            assigning[i] = orbitUsed.get(i);
        }

        ArrayList<int[]> out = new ArrayList<>();
        out.add(partitioning);
        out.add(assigning);
        return out;
    }

    private ArrayList<int[]> swapSubarrays(int[] arr1, int[] arr2, int split){
        int[] out1 = new int[arr1.length];
        int[] out2 = new int[arr2.length];
        for(int i = 0; i < arr1.length; i++){
            if(i < split){
                out1[i] = arr1[i];
                out2[i] = arr2[i];
            }else{
                out2[i] = arr1[i];
                out1[i] = arr2[i];
            }
        }
        ArrayList<int[]> out = new ArrayList<>();
        out.add(out1);
        out.add(out2);
        return out;
    }

    private int[] getIntVariables(String tag, Architecture arch){
        int numVariables = arch.getDecision(tag).getNumberOfVariables();
        int[] variables = new int[numVariables];
        for(int i = 0; i < numVariables; i++){
            variables[i] = Integer.parseInt(arch.getVariable(i).toString());
        }
        return variables;
    }

    private void setIntVariables(String tag, Architecture arch, int[] values){
        ArrayList<Variable> variables = arch.getDecision(tag).getVariables();
        for(int i = 0; i < variables.size(); i++){
            IntegerVariable intVar = (IntegerVariable) variables.get(i);
            intVar.setValue(values[i]);
        }
    }

    @Override
    public int getArity(){ return 2; }
}
