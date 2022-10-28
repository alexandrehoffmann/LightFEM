header_name = "FUNCTION_OPERATORS_HPP"

binaryOperatorAllowed = {}
unaryOperatorAllowed = {}

unaryop_list = ["MINUS", "ABS", "CONJ", "EXP", "LOG", "SQRT", "SIN", "COS", "TAN", "ASIN", "ACOS", "ATAN", "NORM", "TRANSPOSE", "ADJOINT"]
binaryop_list = ["DDOT", "INNER", "OUTER", "CROSS", "SUM", "SUB", "PROD", "DIV"]

binaryop_name = {"SUM":"operator+", "SUB":"operator-", "PROD":"operator*", "DIV":"operator/", "DDOT":"ddot", "INNER":"inner", "OUTER":"outer", "CROSS":"cross"}
unaryop_name = {"MINUS":"operator-", "ABS":"abs", "CONJ":"conj", "EXP":"exp", "LOG":"log", "SQRT":"sqrt", "SIN":"sin", "COS":"cos", "TAN":"tan", "ASIN":"asin", "ACOS":"acos", "ATAN":"atan", "NORM":"norm", "TRANSPOSE":"transpose", "ADJOINT":"adjoint"}

unaryop_for_cpx = {"MINUS":True, "ABS":False, "CONJ":True, "EXP":True, "LOG":True, "SQRT":True, "SIN":True, "COS":True, "TAN":True, "ASIN":True, "ACOS":True, "ATAN":True, "NORM":True, "TRANSPOSE":True, "ADJOINT":True}
unaryop_for_rl = {"MINUS":True, "ABS":True, "CONJ":False, "EXP":True, "LOG":True, "SQRT":True, "SIN":True, "COS":True, "TAN":True, "ASIN":True, "ACOS":True, "ATAN":True, "NORM":True, "TRANSPOSE":True, "ADJOINT":False}

def get_base_name(cpx):
	if cpx:
		return "CpxFunctionExpression"
	return "FunctionExpression"
	
def get_ew_base_name(cpx):
	if cpx:
		return "CpxElementWiseFunctionExpression"
	return "ElementWiseFunctionExpression"

def get_derived_type(derived_expr, tmp_expr):
	if tmp_expr:
		return derived_expr
	return "const " + derived_expr + "&"

def get_derived_slice_type(derived_expr):
	return "typename " + derived_expr + "::ReturnType"

def unary_ret_type(op_name, expr_type, cpx_expr, tmp_expr):
	if cpx_expr:
		return "CpxFunctionUnaryExpression<UnaryOp::" + op_name + ", " + expr_type + ", " + get_derived_type("Expr", tmp_expr) + ">"
	return "FunctionUnaryExpression<UnaryOp::" + op_name + ", " + expr_type + ", " + get_derived_type("Expr", tmp_expr) + ">"
	
def binary_ret_type(op_name, lhs, rhs, cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs, ew_lhs, ew_rhs):
	if ew_lhs:
		if cpx_lhs or cpx_rhs:
			return "CpxElementWiseFunctionBinaryExpression<BinaryOp::" + op_name + ", " + lhs + ", " + get_derived_type("LeftExpr", tmp_lhs) + ", " + rhs + ", " + get_derived_slice_type("RightExpr") + ">"
		return "ElementWiseFunctionBinaryExpression<BinaryOp::" + op_name + ", " + lhs + ", " + get_derived_type("LeftExpr", tmp_lhs) + ", " + rhs + ", " + get_derived_slice_type("RightExpr") + ">"
	if ew_rhs:
		if cpx_lhs or cpx_rhs:
			return "CpxElementWiseFunctionBinaryExpression<BinaryOp::" + op_name + ", " + lhs + ", " + get_derived_slice_type("LeftExpr") + ", " + rhs + ", " + get_derived_type("RightExpr", tmp_rhs) + ">"
		return "ElementWiseFunctionBinaryExpression<BinaryOp::" + op_name + ", " + lhs + ", " + get_derived_slice_type("LeftExpr") + ", " + rhs + ", " + get_derived_type("RightExpr", tmp_rhs) + ">"
	
	if cpx_lhs or cpx_rhs:
		return "CpxFunctionBinaryExpression<BinaryOp::" + op_name + ", " + lhs + ", " + get_derived_type("LeftExpr", tmp_lhs) + ", " + rhs + ", " + get_derived_type("RightExpr", tmp_rhs) + ">"
	return "FunctionBinaryExpression<BinaryOp::" + op_name + ", " + lhs + ", " + get_derived_type("LeftExpr", tmp_lhs) + ", " + rhs + ", " + get_derived_type("RightExpr", tmp_rhs) + ">"

def arg_type(base, derived_expr, cpx_expr, tmp_expr):
	if tmp_expr:
		return get_base_name(cpx_expr) + "<" + base + ", " + derived_expr + ">&&"
	return "const " + get_base_name(cpx_expr) + "<" + base + ", " + derived_expr + ">&"

def ew_arg_type(base, derived_expr, cpx_expr, tmp_expr):
	if tmp_expr:
		return get_ew_base_name(cpx_expr) + "<" + base + ", " + derived_expr + ">&&"
	return "const " + get_ew_base_name(cpx_expr) + "<" + base + ", " + derived_expr + ">&"

def fwd_expr(derived_expr, expr_name, tmp_expr):
	ret = "std::forward<const " + derived_expr + ">( static_cast<const " + derived_expr + "&>(" + expr_name + ") )"
	if tmp_expr:
		ret = "std::forward<" + derived_expr + ">( static_cast<" + derived_expr + "&&>(" + expr_name + ") )"
	return ret
	
def fwd_slice_expr(derived_expr, expr_name, ew_exrp):
	# return "static_cast<" + get_derived_slice_type(derived_expr) + ">(" + expr_name + ")[" + expr_name + ".getMesh()->getElemId(" + ew_exrp + ".getElement())]"
	return expr_name + "[" + expr_name + ".getMesh()->getElemId(" + ew_exrp + ".getElement())]"

with open("operators.hpp", "w") as f:
	f.write("#ifndef " + header_name + "\n")
	f.write("#define " + header_name + "\n")
	f.write("\n")
	f.write("#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>\n")
	f.write("#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>\n")
	f.write("\n")
	f.write("#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>\n")
	f.write("#include <LightFEM/Expression/ElementWise/ElementWiseFunctionBinaryExpression.hpp>\n")
	f.write("#include <LightFEM/Expression/ElementWise/ElementWiseFunctionUnaryExpression.hpp>\n")
	f.write("\n")
	f.write("#include <LightFEM/Expression/Function/FunctionExpression.hpp>\n")
	f.write("#include <LightFEM/Expression/Function/FunctionBinaryExpression.hpp>\n")
	f.write("#include <LightFEM/Expression/Function/FunctionUnaryExpression.hpp>\n")
	f.write("\n")
	f.write("////////////////////////////////////////////////////////////////////////\n")
	f.write("////" +  "min/max operators".center(len("////////////////////////////////////////////////////////////////////////") - 8) + "////\n")
	f.write("////////////////////////////////////////////////////////////////////////\n")
	f.write("\n")
	f.write("template<typename Expr>\n")
	f.write("double min(const FunctionExpression<ExprType::SCALAR, Expr>& expr);\n")
	f.write("\n")
	f.write("template<typename Expr>\n")
	f.write("std::complex< double > min(const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr);\n")
	f.write("\n")
	f.write("template<typename Expr>\n")
	f.write("double max(const FunctionExpression<ExprType::SCALAR, Expr>& expr);\n")
	f.write("\n")
	f.write("template<typename Expr>\n")
	f.write("std::complex< double > max(const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr);\n")
	f.write("\n")
	for op_name in unaryop_list:
		f.write("////////////////////////////////////////////////////////////////////////\n")
		f.write("////" +  (op_name.lower() + " operators").center(len("////////////////////////////////////////////////////////////////////////") - 8) + "////\n")
		f.write("////////////////////////////////////////////////////////////////////////\n")
		f.write("\n")
		for tmp_lhs in [False, True]:
			if unaryop_for_rl[op_name]:
				f.write("template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::" + op_name + ", Type>::isDefined, bool> = true>\n")
				f.write("inline " + unary_ret_type(op_name, "Type", False, tmp_lhs) + " " + unaryop_name[op_name] + "(" + arg_type("Type", "Expr", False, tmp_lhs) + " expr) { return " + unary_ret_type(op_name, "Type", False, tmp_lhs) \
					+ "(" + fwd_expr("Expr", "expr", tmp_lhs) + "); }\n")
				f.write("\n")
			if unaryop_for_cpx[op_name]:
				f.write("template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::" + op_name + ", Type>::isDefined, bool> = true>\n")
				f.write("inline " + unary_ret_type(op_name, "Type", True, tmp_lhs) + " " + unaryop_name[op_name] + "(" + arg_type("Type", "Expr", True, tmp_lhs) + " expr) { return " + unary_ret_type(op_name, "Type", True, tmp_lhs) \
					+ "(" + fwd_expr("Expr", "expr", tmp_lhs) + "); }\n")
				f.write("\n")
	for op_name in binaryop_list:
		f.write("////////////////////////////////////////////////////////////////////////\n")
		f.write("////" +  (op_name.lower() + " operators").center(len("////////////////////////////////////////////////////////////////////////") - 8) + "////\n")
		f.write("////////////////////////////////////////////////////////////////////////\n")
		f.write("\n")
		for cpx_lhs in [False, True]:
			for cpx_rhs in [False, True]:
				for tmp_lhs in [False, True]:
					for tmp_rhs in [False, True]:
						f.write("template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::" + op_name + ", LeftType, RightType>::isDefined, bool> = true>\n")
						f.write("inline " + binary_ret_type(op_name, "LeftType", "RightType", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs, False, False) + " " + binaryop_name[op_name] + "(" + arg_type("LeftType", "LeftExpr", cpx_lhs, tmp_lhs) \
							+ " lhs, " + arg_type("RightType", "RightExpr", cpx_rhs, tmp_rhs) + " rhs) { return " + binary_ret_type(op_name, "LeftType", "RightType", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs, False, False) + "(" \
							+ fwd_expr("LeftExpr", "lhs", tmp_lhs) + ", " + fwd_expr("RightExpr", "rhs", tmp_rhs) + "); }\n")
						f.write("\n")
						f.write("template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::" + op_name + ", LeftType, RightType>::isDefined, bool> = true>\n")
						f.write("inline " + binary_ret_type(op_name, "LeftType", "RightType", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs, True, False) + " " + binaryop_name[op_name] + "(" + ew_arg_type("LeftType", "LeftExpr", cpx_lhs, tmp_lhs) \
							+ " lhs, " + arg_type("RightType", "RightExpr", cpx_rhs, tmp_rhs) + " rhs) { return " + binary_ret_type(op_name, "LeftType", "RightType", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs, True, False) + "(" \
							+ fwd_expr("LeftExpr", "lhs", tmp_lhs) + ", " + fwd_slice_expr("RightExpr", "rhs", "lhs") + "); }\n")
						f.write("\n")
						f.write("template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::" + op_name + ", LeftType, RightType>::isDefined, bool> = true>\n")
						f.write("inline " + binary_ret_type(op_name, "LeftType", "RightType", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs, False, True) + " " + binaryop_name[op_name] + "(" + arg_type("LeftType", "LeftExpr", cpx_lhs, tmp_lhs) \
							+ " lhs, " + ew_arg_type("RightType", "RightExpr", cpx_rhs, tmp_rhs) + " rhs) { return " + binary_ret_type(op_name, "LeftType", "RightType", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs, False, True) + "(" \
							+ fwd_slice_expr("LeftExpr", "lhs", "rhs") + ", " + fwd_expr("RightExpr", "rhs", tmp_rhs) + "); }\n")
						f.write("\n")
	f.write("#include <LightFEM/Expression/Function/operators.tpp>\n")
	f.write("#endif // " + header_name + "\n")
