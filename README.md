## Dynamic Constraint Based Model of Cell Free Protein Synthesis

### Background ###
Cell-free protein expression has emerged as a commonly used application in systems and synthetic biology, and a promising technology for personalized point of care medicine.
Cell-free systems derived from crude cell-extracts have shown remarkable utility as a protein synthesis technology. 
However, if cell-free platforms for on-demand biomanufacturing are to become a reality, the performance limits of these systems must be defined and optimized. Toward this goal, we modeled cell-free protein expression using a dynamic stoichiometric approach constrained to the predictions of previous kinetic modeling. The approach captured the biphasic metabolism and production of a model protein, chloramphenicol acetyltransferase. Flux variability analysis showed that substrate utilization was robust, as CAT production was not affected by the choice of glycolysis or pentose phosphate. Variation of the constraint set revealed central carbon metabolites, specifically upper glycolysis, TCA cycle, and pentose phosphate, to be most effective at training a predictive model; however, fluxes remained unidentifiable. These findings further the knowledge of cell-free protein synthesis and represent a novel tool to inform experimental measurement selection for a variety of metabolic networks, whether in vivo or cell-free.

The dynamic constraint based model is described in the publication:

[Dai D, Horvath N, Vilkhovoy M and J. Varner (2017) Dynamic Constraint Based Model of Cell Free Protein Synthesis]

### Installation
You can download this repository as a zip file, or clone or pull it by using the command:

	git pull https://github.com/varnerlab/publication_cell_free_dFBA_repository.git

or

	git clone https://github.com/varnerlab/publication_cell_free_dFBA_repository.git


Julia must be installed on your machine along with [GLPK](https://github.com/JuliaOpt/GLPK.jl) linear programming solver. Julia can be downloaded/installed on any platform. To install the GLPK program issue the command:

  	julia> Pkg.add("GLPK")

in the Julia REPL.

### Running the model
To run the model, first copy the files from the protein folder you wish to simulate and paste them into the ``Model`` folder. Set the directory to the ``Model`` folder and issue the command ``include("Solve.jl")`` in the Julia REPL.

The ``Solve.jl`` script has several user inputs available:


Output | Description
--- | ---
objective__value | The objective value of the reaction set to be optimized. The default is set to optimize the export of the protein of interest
flux_array | The flux distribution throughout the network. The reaction index can be looked up in Debug.txt
dual_array | Shadow cost
uptake_array | Species array
exit_flag | Status of glpk solver. Solution undefined = 1, solution is feasible = 2, problem has no feasible solution = 4, solution is optimal = 5


#### What each file does

file | description
--- | ---
DataDictionary.jl | Encodes the species and reaction bounds arrays and the objective array. 
Debug.txt | List of reactions and species used to generate the model code.
FluxDriver.jl | Julia interface with the [GLPK](https://github.com/JuliaOpt/GLPK.jl) solver. Users should `NEVER, UNDER ANY CIRCUMSTANCES, EVER` edit this file.
Include.jl | Encodes all the include statements for the project. Should be included at the top of top-level driver scripts.
Network.dat | Stoichiometric array for the model.
<!--Solve.jl | Default top-level driver implementation.-->
Solve_bc.jl | Default top-level driver implementation to solve for the base case.
<!--SolveSingle.jl | Default top-level driver implementation to solve for the single addition sets.-->
Solve_combinations.jl | Default top-level driver implementation to solve for the combination sets.
Solve_svd.jl | Default top-level driver implementation to solve for the singular value decomposition set.
Solve_simulated_annealing.jl | Default top-level driver implementation to solve for the simulated annealing sets.
Bounds.jl | Updates the species and reaction bounds and sets the transcription and translation rates.
TXTLDictionary.jl | Encodes the cell-free protein synthesis parameters. Data is stored in a [Julia dictionary](http://docs.julialang.org/en/stable/stdlib/collections/?highlight=dict#Base.Dict) type and can be accessed through the appropriate key.
calculate_constraints.jl | Calculates species constraints, the maximum rate of accumulation and depletion, based on concentrations, measurement data, and user definition. 



### Model code and parameter ensemble
The constraint based wwas implemented in [Julia](http://julialang.org) and solved using the GLPK routine of the [GLPK Package](https://github.com/JuliaOpt/GLPK.jl). The model code is freely available under an [MIT software license](https://opensource.org/licenses/MIT).

To run the base case model, run ``Solve_bc.jl`` The model stoichiometric matrix is encoded in ``Network.dat`` which is called by the ``DataDictionary.jl`` to formulate the optimization problem along with the constraint initialization. The transcription and translation constrainst are specified in ``TXTLDictionary.jl`` which is called by ``Bounds.jl`` that incorporated all the flux constraints. Metabolite constraints are incorporated by ``calculate_constraints.jl`` which requires metabolite concentration. Finally, the optimization problem is solved with ``FluxDriver.jl`` that set up the optimization problem to be solved with GLPK. To execute this script, issue the command:

``julia> include("Solve_bc.jl")``

__Prerequisites__: [Julia](http://julialang.org) and the [GLPK Package](https://github.com/JuliaOpt/GLPK.jl) must be installed on your computer before the model equations can be solved.


