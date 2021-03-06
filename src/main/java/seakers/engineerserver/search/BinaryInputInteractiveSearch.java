/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seakers.engineerserver.search;


import com.google.gson.*;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.util.TypedProperties;
import seakers.engineerserver.thriftinterface.BinaryInputArchitecture;

import java.util.ArrayList;

/**
 *
 * @author antoni
 */
public class BinaryInputInteractiveSearch extends AbstractInteractiveSearch {

    public BinaryInputInteractiveSearch(Algorithm alg, TypedProperties properties, String id) {
        super(alg, properties, id);
    }

    @Override
    public JsonElement getJSONArchitecture(Solution architecture) {
        BinaryInputArchitecture jsonArch = new BinaryInputArchitecture();
        jsonArch.inputs = new ArrayList<>();
        jsonArch.outputs = new ArrayList<>();
        for (int j = 1; j < architecture.getNumberOfVariables(); ++j) {
            BinaryVariable var = (BinaryVariable)architecture.getVariable(j);
            boolean binaryVal = var.get(0);
            jsonArch.inputs.add(binaryVal);
        }
        jsonArch.outputs.add(-architecture.getObjective(0));
        jsonArch.outputs.add(architecture.getObjective(1));

        Gson gson = new GsonBuilder().create();
        return gson.toJsonTree(jsonArch);
    }

}
