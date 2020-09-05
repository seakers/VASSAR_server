/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package seakers.engineerserver.search.problems.Scheduling;

import org.moeaframework.core.Solution;
import seakers.architecture.pattern.ArchitecturalDecision;
import seakers.architecture.pattern.Permuting;
import seakers.architecture.Architecture;

import java.util.ArrayList;

/**
 * This class creates a solution for the problem consisting of an assigning
 * pattern of instruments to orbits and a combining pattern for the number of
 * satellites per orbit. Assigning instruments from the left hand side to orbits
 * on the right hand side
 *
 * @author nozomi
 */
public class SchedulingArchitecture extends Architecture {

    private static final long serialVersionUID = 8776271523867355732L;

    /**
     * Tag used for the permuting decision
     */
    private static final String permuteTag = "perm";

    private boolean alreadyEvaluated;

    //Constructors
    /**
     * Creates an empty architecture with a default number of satellites.
     * Default value is the first value in the array given as
     * alternativesForNumberOfSatellites
     *
     * @param numberOfMissions
     * @param numberOfObjectives
     */
    public SchedulingArchitecture(int numberOfMissions, int numberOfObjectives) {
        super(numberOfObjectives, 0,
                createDecisions(numberOfMissions));
        this.alreadyEvaluated = false;
    }

    private static ArrayList<ArchitecturalDecision> createDecisions(int numberOfMissions) {
        ArrayList<ArchitecturalDecision> out = new ArrayList<>();
        out.add(new Permuting(numberOfMissions, permuteTag));
        return out;
    }

    /**
     * makes a copy solution from the input solution
     *
     * @param solution
     */
    private SchedulingArchitecture(Solution solution) {
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
        return new SchedulingArchitecture(this);
    }
}
