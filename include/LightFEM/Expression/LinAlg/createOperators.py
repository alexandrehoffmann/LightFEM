header_name = "LINALG_OPERATORS_HPP"

binaryOperatorAllowed = {}
unaryOperatorAllowed = {}

unaryop_list = ["MINUS", "ABS", "CONJ", "EXP", "LOG", "SQRT", "SIN", "COS", "TAN", "ASIN", "ACOS", "ATAN", "NORM", "TRANSPOSE", "ADJOINT"]
binaryop_list = ["DDOT", "INNER", "OUTER", "CROSS", "SUM", "SUB", "PROD", "DIV"]

type_list = ["RK4_TENSOR", "RK2_TENSOR", "MATRIX", "VECTOR", "SCALAR"]
base_list = {"RK4_TENSOR":"RankFourTensorExpression", "RK2_TENSOR":"RankTwoTensorExpression", "MATRIX":"MatrixExpression", "VECTOR":"VectorExpression", "SCALAR":"ScalarExpression"}

binaryop_name = {"SUM":"operator+", "SUB":"operator-", "PROD":"operator*", "DIV":"operator/", "DDOT":"ddot", "INNER":"inner", "OUTER":"outer", "CROSS":"cross"}
unaryop_name = {"MINUS":"operator-", "ABS":"abs", "CONJ":"conj", "EXP":"exp", "LOG":"log", "SQRT":"sqrt", "SIN":"sin", "COS":"cos", "TAN":"tan", "ASIN":"asin", "ACOS":"acos", "ATAN":"atan", "NORM":"norm", "TRANSPOSE":"transpose", "ADJOINT":"adjoint"}

def get_base_name(base, cpx):
	if cpx:
		return "Cpx" + base_list[base]
	return base_list[base]

def get_derived_type(derived_expr, tmp_expr):
	if tmp_expr:
		return derived_expr
	return "const " + derived_expr + "&"

def unary_ret_type(op_name, expr_type, cpx_expr, tmp_expr):
	if cpx_expr:
		return "CpxUnaryExpression<UnaryOp::" + op_name + ", ExprType::" + expr_type + ", " + get_derived_type("Expr", tmp_expr) + ">"
	return "UnaryExpression<UnaryOp::" + op_name + ", ExprType::" + expr_type + ", " + get_derived_type("Expr", tmp_expr) + ">"
	
def binary_ret_type(op_name, lhs, rhs, cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs):
	if cpx_lhs or cpx_rhs:
		return "CpxBinaryExpression<BinaryOp::" + op_name + ", ExprType::" + lhs + ", " + get_derived_type("LeftExpr", tmp_lhs) + ", ExprType::" + rhs + ", " + get_derived_type("RightExpr", tmp_rhs) + ">"
	return "BinaryExpression<BinaryOp::" + op_name + ", ExprType::" + lhs + ", " + get_derived_type("LeftExpr", tmp_lhs) + ", ExprType::" + rhs + ", " + get_derived_type("RightExpr", tmp_rhs) + ">"

def arg_type(base, derived_expr, cpx_expr, tmp_expr):
	if tmp_expr:
		return get_base_name(base, cpx_expr) + "<" + derived_expr + ">&&"
	return "const " + get_base_name(base, cpx_expr) + "<" + derived_expr + ">&"
	
def fwd_expr(derived_expr, expr_name, tmp_expr):
	ret = "std::forward<const " + derived_expr + ">( static_cast<const " + derived_expr + "&>(" + expr_name + ") )"
	if tmp_expr:
		ret = "std::forward<" + derived_expr + ">( static_cast<" + derived_expr + "&&>(" + expr_name + ") )"
	return ret

for unaryop in unaryop_list:
	for expr in type_list:
		unaryOperatorAllowed[(unaryop, expr)] = (False, False)
for binaryop in binaryop_list:
	for lhs in type_list:
		for rhs in type_list:
			binaryOperatorAllowed[(binaryop, lhs, rhs)] = False

unaryOperatorAllowed[("MINUS","RK4_TENSOR")] = (True, True)
unaryOperatorAllowed[("MINUS","RK2_TENSOR")] = (True, True)
unaryOperatorAllowed[("MINUS","MATRIX")] = (True, True)
unaryOperatorAllowed[("MINUS","VECTOR")] = (True, True)
unaryOperatorAllowed[("MINUS","SCALAR")] = (True, True)

unaryOperatorAllowed[("ABS","SCALAR")] = (True, False)

unaryOperatorAllowed[("CONJ","RK4_TENSOR")] = (False, True)
unaryOperatorAllowed[("CONJ","RK2_TENSOR")] = (False, True)
unaryOperatorAllowed[("CONJ","MATRIX")] = (False, True)
unaryOperatorAllowed[("CONJ","VECTOR")] = (False, True)
unaryOperatorAllowed[("CONJ","SCALAR")] = (False, True)

unaryOperatorAllowed[("EXP","SCALAR")] = (True, True)
unaryOperatorAllowed[("LOG","SCALAR")] = (True, True)
unaryOperatorAllowed[("SQRT","SCALAR")] = (True, True)
unaryOperatorAllowed[("SIN","SCALAR")] = (True, True)
unaryOperatorAllowed[("COS","SCALAR")] = (True, True)
unaryOperatorAllowed[("TAN","SCALAR")] = (True, True)
unaryOperatorAllowed[("ASIN","SCALAR")] = (True, True)
unaryOperatorAllowed[("ACOS","SCALAR")] = (True, True)
unaryOperatorAllowed[("ATAN","SCALAR")] = (True, True)

unaryOperatorAllowed[("NORM","RK2_TENSOR")] = (True, True)
unaryOperatorAllowed[("NORM","MATRIX")] = (True, True)
unaryOperatorAllowed[("NORM","VECTOR")] = (True, True)

unaryOperatorAllowed[("TRANSPOSE","RK4_TENSOR")] = (True, True)
unaryOperatorAllowed[("TRANSPOSE","RK2_TENSOR")] = (True, True)
unaryOperatorAllowed[("TRANSPOSE","MATRIX")] = (True, True)

unaryOperatorAllowed[("ADJOINT","RK4_TENSOR")] = (False, True)
unaryOperatorAllowed[("ADJOINT","RK2_TENSOR")] = (False, True)
unaryOperatorAllowed[("ADJOINT","MATRIX")] = (False, True)
		
binaryOperatorAllowed[("SUM", "RK4_TENSOR", "RK4_TENSOR")] = True
binaryOperatorAllowed[("SUM", "RK2_TENSOR", "RK2_TENSOR")] = True
binaryOperatorAllowed[("SUM", "MATRIX", "MATRIX")] = True
binaryOperatorAllowed[("SUM", "VECTOR", "VECTOR")] = True
binaryOperatorAllowed[("SUM", "SCALAR", "SCALAR")] = True

binaryOperatorAllowed[("SUB", "RK4_TENSOR", "RK4_TENSOR")] = True
binaryOperatorAllowed[("SUB", "RK2_TENSOR", "RK2_TENSOR")] = True
binaryOperatorAllowed[("SUB", "MATRIX", "MATRIX")] = True
binaryOperatorAllowed[("SUB", "VECTOR", "VECTOR")] = True
binaryOperatorAllowed[("SUB", "SCALAR", "SCALAR")] = True

binaryOperatorAllowed[("PROD", "RK4_TENSOR", "SCALAR")] = True
binaryOperatorAllowed[("PROD", "RK2_TENSOR", "SCALAR")] = True
binaryOperatorAllowed[("PROD", "MATRIX", "MATRIX")] = True
binaryOperatorAllowed[("PROD", "MATRIX", "VECTOR")] = True
binaryOperatorAllowed[("PROD", "MATRIX", "SCALAR")] = True
binaryOperatorAllowed[("PROD", "VECTOR", "SCALAR")] = True
binaryOperatorAllowed[("PROD", "SCALAR", "RK4_TENSOR")] = True
binaryOperatorAllowed[("PROD", "SCALAR", "RK2_TENSOR")] = True
binaryOperatorAllowed[("PROD", "SCALAR", "MATRIX")] = True
binaryOperatorAllowed[("PROD", "SCALAR", "VECTOR")] = True
binaryOperatorAllowed[("PROD", "SCALAR", "SCALAR")] = True

binaryOperatorAllowed[("DIV", "RK4_TENSOR", "SCALAR")] = True
binaryOperatorAllowed[("DIV", "RK2_TENSOR", "SCALAR")] = True
binaryOperatorAllowed[("DIV", "MATRIX", "SCALAR")] = True
binaryOperatorAllowed[("DIV", "VECTOR", "SCALAR")] = True
binaryOperatorAllowed[("DIV", "SCALAR", "SCALAR")] = True

binaryOperatorAllowed[("DDOT", "RK4_TENSOR", "RK2_TENSOR")] = True
binaryOperatorAllowed[("DDOT", "RK4_TENSOR", "MATRIX")] = True

binaryOperatorAllowed[("INNER", "RK2_TENSOR", "RK2_TENSOR")] = True
binaryOperatorAllowed[("INNER", "MATRIX", "MATRIX")] = True
binaryOperatorAllowed[("INNER", "VECTOR", "VECTOR")] = True

binaryOperatorAllowed[("OUTER", "VECTOR", "VECTOR")] = True

binaryOperatorAllowed[("CROSS", "VECTOR", "VECTOR")] = True

with open("operators.hpp", "w") as f:
	f.write("#ifndef " + header_name + "\n")
	f.write("#define " + header_name + "\n")
	f.write("\n")
	f.write("#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>\n")
	f.write("#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>\n")
	f.write("\n")
	for lhs in type_list:
		f.write("////////////////////////////////////////////////////////////////////////\n")
		f.write("////" + (lhs.lower() + " operators").center(len("////////////////////////////////////////////////////////////////////////") - 8) + "////\n")
		f.write("////////////////////////////////////////////////////////////////////////\n")
		f.write("\n")
		for op_name in unaryop_list:
			if unaryOperatorAllowed[(op_name, lhs)][0]:
				for tmp_lhs in [False, True]:
					f.write("template<typename Expr>\n")
					f.write("inline " + unary_ret_type(op_name, lhs, False, tmp_lhs) + " " + unaryop_name[op_name] + "(" + arg_type(lhs, "Expr", False, tmp_lhs) + " expr) { return " + unary_ret_type(op_name, lhs, False, tmp_lhs) \
						+ "(" + fwd_expr("Expr", "expr", tmp_lhs) + "); }\n")
					f.write("\n")
			if unaryOperatorAllowed[(op_name, lhs)][1]:
				for tmp_lhs in [False, True]:
					f.write("template<typename Expr>\n")
					f.write("inline " + unary_ret_type(op_name, lhs, True, tmp_lhs) + " " + unaryop_name[op_name] + "(" + arg_type(lhs, "Expr", True, tmp_lhs) + " expr) { return " + unary_ret_type(op_name, lhs, True, tmp_lhs) \
						+ "(" + fwd_expr("Expr", "expr", tmp_lhs) + "); }\n")
					f.write("\n")
		for op_name in binaryop_list:
			for rhs in type_list:
				if binaryOperatorAllowed[(op_name, lhs, rhs)]:
					for cpx_lhs in [False, True]:
						for cpx_rhs in [False, True]:
							for tmp_lhs in [False, True]:
								for tmp_rhs in [False, True]:
									f.write("template<typename LeftExpr, typename RightExpr>\n")
									f.write("inline " + binary_ret_type(op_name, lhs, rhs, cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs) + " " + binaryop_name[op_name] + "(" + arg_type(lhs, "LeftExpr", cpx_lhs, tmp_lhs) + " lhs, " \
										+ arg_type(rhs, "RightExpr", cpx_rhs, tmp_rhs) + " rhs) { return " + binary_ret_type(op_name, lhs, rhs, cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs) + "(" + fwd_expr("LeftExpr", "lhs", tmp_lhs) + ", " \
										+ fwd_expr("RightExpr", "rhs", tmp_rhs) + "); }\n")
									f.write("\n")
	f.write("#endif // " + header_name + "\n")
