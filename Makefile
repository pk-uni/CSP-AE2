idf : 
	minizinc src/cp/infectious-defence.mzn src/cp/data/instance-$(DATA).dzn

base:
	minizinc src/cp/base/integer.mzn src/cp/data/instance-$(DATA).dzn