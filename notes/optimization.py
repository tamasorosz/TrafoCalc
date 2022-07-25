import json

from math import inf
from artap.algorithm_genetic import NSGAII
from artap.problem import Problem

from importlib_resources import files
from src.two_winding_model import TransformerDesign, TwoWindingModel
from src.models import IndependentVariables


class TransformerOptimizationProblem(Problem):
    def set(self):
        self.name = "Optimization of the 10 MVA validation problem in the case of two different scenarios."

        self.parameters = [
            {"name": "rc", "bounds": [180.0, 250.0]},
            {"name": "bc", "bounds": [1.5, 1.7]},
            {"name": "j_in", "bounds": [2.0, 3.0]},
            {"name": "j_ou", "bounds": [2.0, 3.0]},
            {"name": "h_in", "bounds": [800., 1400.]},
            {"name": "m_gap", "bounds": [20.0, 60.0]},
        ]

        self.costs = [{"name": "TOC", "criteria": "minimize"}]

        path = files("data").joinpath("10MVA_example.json")

        with open(path) as json_file:
            self.transformer_data = json.load(json_file)

    def individual_status(self, x):
        """
        Prints out the selected optimization parameters of the current individual.
        :param x: the list of the optimization parameters
        """
        x1 = [round(xi, 2) for xi in x]
        print("called with", x1, end=" ")
        assert len(x) == 6

    def evaluate(self, individual):
        x = individual.vector
        self.individual_status(x)

        try:
            transformer = TransformerDesign.from_dict(self.transformer_data)
            transformer.design_params = IndependentVariables(rc=x[0], bc=x[1], j_in=x[2], j_ou=x[3], h_in=x[4],
                                                             m_gap=x[5])

            trafo_model = TwoWindingModel(input=transformer)
            trafo_model.calculate(is_sc=False)

            # FEM calculation
            trafo_model.fem_simulation(detailed_output=False)

            if not transformer.required.check_sci_requrements(trafo_model.results.fem_based_sci):
                return [inf]

            # f1 and the modified f2 and f3 measures needs only one evaluation
            res = trafo_model.results.capitalized_cost

            print("Capitalized cost:", round(res, 0))
            print("DONE")

            return [res]

        except:
            print("FAILED")
            return [inf]


if __name__ == "__main__":

    # Perform the optimization iterating over 100 times on 100 individuals.
    problem = TransformerOptimizationProblem()
    algorithm = NSGAII(problem)
    algorithm.options["max_population_number"] = 30
    algorithm.options["max_population_size"] = 30
    try:
        algorithm.run()
        res = problem.individuals[-1]
        print("OPTIMIZATION RESULT:")
        print(res.vector)
        print(res.costs)
    except KeyboardInterrupt:
        pass
