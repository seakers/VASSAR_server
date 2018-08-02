/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package search.problems.PartitioningAndAssigning;

import org.moeaframework.core.Solution;
import seak.architecture.pattern.ArchitecturalDecision;
import seak.architecture.pattern.Assigning;
import seak.architecture.pattern.Partitioning;

import java.util.ArrayList;


public class PartitioningAndAssigningArchitecture extends seak.architecture.Architecture {

    private static final long serialVersionUID = 8776271523867355732L;

    /**
     * Tag used for the assigning decision
     */
    private static final String assignTag = "orbitAssignment";

    /**
     * Tag used for the combining decision
     */
    private static final String partitionTag = "instrumentPartitioning";

    private boolean alreadyEvaluated;

    //Constructors
    /**
     * Creates an empty architecture with a default number of satellites.
     * Default value is the first value in the array given as
     * alternativesForNumberOfSatellites
     *
     * @param numberOfInstruments
     * @param numberOfOrbits
     * @param numberOfObjectives
     */
    public PartitioningAndAssigningArchitecture(int numberOfInstruments, int numberOfOrbits, int numberOfObjectives) {
        super(numberOfObjectives, 1, createDecisions(numberOfInstruments, numberOfOrbits));
        this.alreadyEvaluated = false;
    }

    private static ArrayList<ArchitecturalDecision> createDecisions(int numberOfInstruments, int numberOfOrbits) {
        ArrayList<ArchitecturalDecision> out = new ArrayList<>();
        out.add(new Partitioning(numberOfInstruments, partitionTag));
        out.add(new Assigning(numberOfInstruments, numberOfOrbits, assignTag));
        return out;
    }

    /**
     * makes a copy solution from the input solution
     *
     * @param solution
     */
    private PartitioningAndAssigningArchitecture(Solution solution) {
        super(solution);
    }

    public void setAlreadyEvaluated(boolean alreadyEvaluated) {
        this.alreadyEvaluated = alreadyEvaluated;
    }

    public boolean getAlreadyEvaluated() {
        return this.alreadyEvaluated;
    }

    @Override
    public Solution copy() {
        return new PartitioningAndAssigningArchitecture(this);
    }
}
