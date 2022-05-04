from importlib_resources import files
from src.two_winding_model import TransformerDesign, TwoWindingModel

path = files("data").joinpath("1250kVA_sc_transformer.json")

import json

with open(path) as json_file:
    data = json.load(json_file)

transformer = TransformerDesign.from_dict(data)

trafo_model = TwoWindingModel(input=transformer)
trafo_model.calculate(is_sc=True)

# FEM calculation
trafo_model.fem_simulation()

print(trafo_model.hv_winding)
print(trafo_model.lv_winding)
print(trafo_model.results)
print('analytical sci:', trafo_model.results.sci)