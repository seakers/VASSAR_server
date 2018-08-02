package search.problems.PartitioningAndAssigning.operators;

import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.Variation;
import rbsa.eoss.local.BaseParams;
import seak.architecture.Architecture;
import seak.architecture.operators.IntegerUM;
import seak.architecture.util.IntegerVariable;
import search.problems.PartitioningAndAssigning.PartitioningAndAssigningArchitecture;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Collection;


public class PartitioningAndAssigningMutation extends IntegerUM{

    protected double probability;
    protected BaseParams params;

    public PartitioningAndAssigningMutation(double probability, BaseParams params){
        super(probability);
        this.probability = probability;
        this.params = params;
    }

    @Override
    public Solution[] evolve(Solution[] parents) {
        Solution[] out = super.evolve(parents);
        out[0] = repair(out[0]);
        return out;
    }

    private Solution repair(Solution solution){

        Architecture arch = (PartitioningAndAssigningArchitecture) solution;
        int[] partitioning = getIntVariables("instrumentPartitioning", arch);
        int[] assigning = getIntVariables("orbitAssignment", arch);

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
        }

        setIntVariables("instrumentPartitioning", newArch, newPartitioning);
        setIntVariables("orbitAssignment", newArch, newAssigning);
        return newArch;
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
}
