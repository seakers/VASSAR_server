


/**
 *  The available types in Thrift are:
 *
 *  bool        Boolean, one byte
 *  i8 (byte)   Signed 8-bit integer
 *  i16         Signed 16-bit integer
 *  i32         Signed 32-bit integer
 *  i64         Signed 64-bit integer
 *  double      64-bit floating point value
 *  string      String
 *  binary      Blob (byte array)
 *  map<t1,t2>  Map from one type to another
 *  list<t1>    Ordered list of one type
 *  set<t1>     Set of unique elements of one type
 *
 */


namespace java javaInterface
namespace py pyInterface

typedef i32 int




/**
 * Structs are the basic complex data structures. They are comprised of fields
 * which each have an integer identifier, a type, a symbolic name, and an
 * optional default value.
 *
 * Fields can be declared "optional", which ensures they will not be included
 * in the serialized output if they aren't set.  Note that this requires some
 * manual management in some languages.
 */
 
 
 
struct BinaryInputArchitecture {
    1: int id,
    2: list<bool> inputs,
    3: list<double> outputs
}

struct ObjectiveSatisfaction {
    1: string objective_name,
    2: double satisfaction,
    3: double weight
}

struct SubscoreInformation {
    1: string name,
    2: string description,
    3: double value,
    4: double weight,
    5: list<SubscoreInformation> subscores
}

struct MissionCostInformation {
    1: string orbit_name,
    2: list<string> payload,
    3: string launch_vehicle,
    4: double total_mass,
    5: double total_power,
    6: double total_cost,
    7: map<string, double> mass_budget,
    8: map<string, double> power_budget,
    9: map<string, double> cost_budget
}

struct SubobjectiveDetails {
    1: string param,
    2: list<string> attr_names,
    3: list<list<string>> attr_values,
    4: list<double> scores,
    5: list<string> taken_by,
    6: list<list<string>> justifications
}


service VASSARInterface {

    /**
    * A method definition looks like C code. It has a return type, arguments,
    * and optionally a list of exceptions that it may throw. Note that argument
    * lists and exception lists are specified using the exact same syntax as
    * field lists in struct or exception definitions.
    */

    void ping(),

    BinaryInputArchitecture eval(1:list<bool> inputs),

    list<BinaryInputArchitecture> runLocalSearch(1:list<bool> inputs),

    list<string> getOrbitList(),

    list<string> getInstrumentList(),

    list<string> getObjectiveList(),

    list<string> getInstrumentsForObjective(1:string objective),

    list<string> getInstrumentsForPanel(1:string panel),

    list<string> getCritique(1:list<bool> inputs),

    list<ObjectiveSatisfaction> getArchitectureScoreExplanation(1:list<bool> arch),

    list<ObjectiveSatisfaction> getPanelScoreExplanation(1:list<bool> arch, 2:string panel),

    list<ObjectiveSatisfaction> getObjectiveScoreExplanation(1:list<bool> arch, 2:string objective),

    oneway void startGA(1:list<BinaryInputArchitecture> dataset, 2:string username),

    list<SubscoreInformation> getArchScienceInformation(1:BinaryInputArchitecture arch),

    list<MissionCostInformation> getArchCostInformation(1:BinaryInputArchitecture arch),

    SubobjectiveDetails getSubscoreDetails(1:BinaryInputArchitecture arch, 2:string subobj)
}


