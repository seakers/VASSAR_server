/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seakers.engineerserver.search;


import com.google.gson.*;
import seakers.engineerserver.thriftinterface.DiscreteInputArchitecture;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Solution;
import org.moeaframework.util.TypedProperties;
import seakers.architecture.util.IntegerVariable;

import java.util.ArrayList;

/**
 *
 * @author antoni
 */
public class DiscreteInputInteractiveSearch extends AbstractInteractiveSearch {

    public DiscreteInputInteractiveSearch(Algorithm alg, TypedProperties properties, String username) {
        super(alg, properties, username);
    }

    @Override
    public JsonElement getJSONArchitecture(Solution architecture) {
        DiscreteInputArchitecture json_arch = new DiscreteInputArchitecture();
        json_arch.inputs = new ArrayList<>();
        json_arch.outputs = new ArrayList<>();
        for (int j = 0; j < architecture.getNumberOfVariables(); ++j) {
            IntegerVariable var = (IntegerVariable)architecture.getVariable(j);
            int val = var.getValue();
            json_arch.inputs.add(val);
        }
        json_arch.outputs.add(-architecture.getObjective(0));
        json_arch.outputs.add(architecture.getObjective(1));

        Gson gson = new GsonBuilder().create();
        return gson.toJsonTree(json_arch);
    }

}
