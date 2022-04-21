header_name = "FINITE_ELEMENT_FUNCTION_OPERATORS_HPP"

binaryOperatorAllowed = {}
unaryOperatorAllowed = {}

unaryop_list = ["CONJ", "TRANSPOSE", "ADJOINT", "MINUS"]
binaryop_list = ["SUM", "SUB"]

binaryop_name = {"SUM":"operator+", "SUB":"operator-", "PROD":"operator*", "DIV":"operator/", "DDOT":"ddot", "INNER":"inner", "OUTER":"outer", "CROSS":"cross"}
unaryop_name = {"MINUS":"operator-", "CONJ":"conj", "NORM":"norm", "TRANSPOSE":"transpose", "ADJOINT":"adjoint"}

unaryop_for_cpx = {"MINUS":True, "CONJ":True, "NORM":True, "TRANSPOSE":True, "ADJOINT":True}
unaryop_for_rl = {"MINUS":True, "CONJ":False, "NORM":True, "TRANSPOSE":True, "ADJOINT":False}

def get_base_name(cpx):
	if cpx:
		return "CpxFiniteElementFunctionExpression"
	return "FiniteElementFunctionExpression"

def get_derived_type(derived_expr, tmp_expr):
	if tmp_expr:
		return derived_expr
	return "const " + derived_expr + "&"

def unary_ret_type(op_name, expr_type, cpx_expr, tmp_expr):
	if cpx_expr:
		return "CpxFiniteElementFunctionUnaryExpression<UnaryOp::" + op_name + ", " + expr_type + ", " + get_derived_type("Expr", tmp_expr) + ">"
	return "FiniteElementFunctionUnaryExpression<UnaryOp::" + op_name + ", " + expr_type + ", " + get_derived_type("Expr", tmp_expr) + ">"
	
def binary_ret_type(op_name, lhs, rhs, cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs):
	if cpx_lhs or cpx_rhs:
		return "CpxFiniteElementFunctionBinaryExpression<BinaryOp::" + op_name + ", " + lhs + ", " + get_derived_type("LeftExpr", tmp_lhs) + ", " + rhs + ", " + get_derived_type("RightExpr", tmp_rhs) + ">"
	return "FiniteElementFunctionBinaryExpression<BinaryOp::" + op_name + ", " + lhs + ", " + get_derived_type("LeftExpr", tmp_lhs) + ", " + rhs + ", " + get_derived_type("RightExpr", tmp_rhs) + ">"
	
def binary_ret_type_lhscst(op_name, lhs, rhs, cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs):
	dtype = "Const"
	if cpx_lhs:
		dtype = "Cpx" + dtype
	if cpx_lhs or cpx_rhs:
		return "CpxFiniteElementFunctionBinaryExpression<BinaryOp::" + op_name + ", " + lhs + ", " + get_derived_type(dtype, tmp_lhs) + ", " + rhs + ", " + get_derived_type("RightExpr", tmp_rhs) + ">"
	return "FiniteElementFunctionBinaryExpression<BinaryOp::" + op_name + ", " + lhs + ", " + get_derived_type(dtype, tmp_lhs) + ", " + rhs + ", " + get_derived_type("RightExpr", tmp_rhs) + ">"

def binary_ret_type_rhscst(op_name, lhs, rhs, cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs):
	dtype = "Const"
	if cpx_rhs:
		dtype = "Cpx" + dtype
	if cpx_lhs or cpx_rhs:
		return "CpxFiniteElementFunctionBinaryExpression<BinaryOp::" + op_name + ", " + lhs + ", " + get_derived_type("LeftExpr", tmp_lhs) + ", " + rhs + ", " + get_derived_type(dtype, tmp_rhs) + ">"
	return "FiniteElementFunctionBinaryExpression<BinaryOp::" + op_name + ", " + lhs + ", " + get_derived_type("LeftExpr", tmp_lhs) + ", " + rhs + ", " + get_derived_type(dtype, tmp_rhs) + ">"

def arg_type(base, derived_expr, cpx_expr, tmp_expr):
	if tmp_expr:
		return get_base_name(cpx_expr) + "<" + base + ", " + derived_expr + ">&&"
	return "const " + get_base_name(cpx_expr) + "<" + base + ", " + derived_expr + ">&"

def const_type(cpx_expr, tmp_expr):
	type_expr = "Const"
	if cpx_expr:
		type_expr = "Cpx" + type_expr
	if tmp_expr:
		return type_expr + "&&"
	return "const " + type_expr + "&"

def fwd_expr(derived_expr, expr_name, tmp_expr):
	ret = "std::forward<const " + derived_expr + ">( static_cast<const " + derived_expr + "&>(" + expr_name + ") )"
	if tmp_expr:
		ret = "std::forward<" + derived_expr + ">( static_cast<" + derived_expr + "&&>(" + expr_name + ") )"
	return ret

with open("operators.hpp", "w") as f:
	f.write("#ifndef " + header_name + "\n")
	f.write("#define " + header_name + "\n")
	f.write("\n")
	f.write("#include <Expression/LinAlg/BinaryExpression.hpp>\n")
	f.write("#include <Expression/LinAlg/UnaryExpression.hpp>\n")
	f.write("\n")
	f.write("#include <Expression/FiniteElementFunction/FiniteElementFunctionExpression.hpp>\n")
	f.write("#include <Expression/FiniteElementFunction/FiniteElementFunctionBinaryExpression.hpp>\n")
	f.write("#include <Expression/FiniteElementFunction/FiniteElementFunctionUnaryExpression.hpp>\n")
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
						f.write("inline " + binary_ret_type(op_name, "LeftType", "RightType", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs) + " " + binaryop_name[op_name] + "(" + arg_type("LeftType", "LeftExpr", cpx_lhs, tmp_lhs) \
							+ " lhs, " + arg_type("RightType", "RightExpr", cpx_rhs, tmp_rhs) + " rhs) { return " + binary_ret_type(op_name, "LeftType", "RightType", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs) + "(" \
							+ fwd_expr("LeftExpr", "lhs", tmp_lhs) + ", " + fwd_expr("RightExpr", "rhs", tmp_rhs) + "); }\n")
						f.write("\n")
	f.write("////////////////////////////////////////////////////////////////////////\n")
	f.write("////                         prod operators                         ////\n")
	f.write("////////////////////////////////////////////////////////////////////////\n")
	f.write("\n")
	op_name = "PROD"
	for cpx_lhs in [False, True]:
		for cpx_rhs in [False, True]:
			for tmp_lhs in [False, True]:
				for tmp_rhs in [False, True]:
					const_name  = "Const"
					if cpx_lhs:
						const_name = "Cpx" + const_name
					f.write("template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::" + op_name + ", ExprType::SCALAR, RightType>::isDefined, bool> = true>\n")
					f.write("inline " + binary_ret_type_lhscst(op_name, "ExprType::SCALAR", "RightType", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs) + " " + binaryop_name[op_name] + "(" + const_type(cpx_lhs, tmp_lhs) \
						+ " lhs, " + arg_type("RightType", "RightExpr", cpx_rhs, tmp_rhs) + " rhs) { return " + binary_ret_type_lhscst(op_name, "ExprType::SCALAR", "RightType", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs) + "(" \
						+ fwd_expr(const_name, "lhs", tmp_lhs) + ", " + fwd_expr("RightExpr", "rhs", tmp_rhs) + "); }\n")
					f.write("\n")
	for cpx_lhs in [False, True]:
		for cpx_rhs in [False, True]:
			for tmp_lhs in [False, True]:
				for tmp_rhs in [False, True]:
					const_name  = "Const"
					if cpx_rhs:
						const_name = "Cpx" + const_name
					f.write("template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::" + op_name + ", LeftType, ExprType::SCALAR>::isDefined, bool> = true>\n")
					f.write("inline " + binary_ret_type_rhscst(op_name, "LeftType", "ExprType::SCALAR", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs) + " " + binaryop_name[op_name] + "(" + arg_type("LeftType", "LeftExpr", cpx_lhs, tmp_lhs) \
						+ " lhs, " + const_type(cpx_rhs, tmp_rhs) + " rhs) { return " + binary_ret_type_rhscst(op_name, "LeftType", "ExprType::SCALAR", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs) + "(" \
						+ fwd_expr("LeftExpr", "lhs", tmp_lhs) + ", " + fwd_expr(const_name, "rhs", tmp_rhs) + "); }\n")
					f.write("\n")
	op_name = "DIV"
	for cpx_lhs in [False, True]:
		for cpx_rhs in [False, True]:
			for tmp_lhs in [False, True]:
				for tmp_rhs in [False, True]:
					const_name  = "Const"
					if cpx_rhs:
						const_name = "Cpx" + const_name
					f.write("template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::" + op_name + ", LeftType, ExprType::SCALAR>::isDefined, bool> = true>\n")
					f.write("inline " + binary_ret_type_rhscst(op_name, "LeftType", "ExprType::SCALAR", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs) + " " + binaryop_name[op_name] + "(" + arg_type("LeftType", "LeftExpr", cpx_lhs, tmp_lhs) \
						+ " lhs, " + const_type(cpx_rhs, tmp_rhs) + " rhs) { return " + binary_ret_type_rhscst(op_name, "LeftType", "ExprType::SCALAR", cpx_lhs, cpx_rhs, tmp_lhs, tmp_rhs) + "(" \
						+ fwd_expr("LeftExpr", "lhs", tmp_lhs) + ", " + fwd_expr(const_name, "rhs", tmp_rhs) + "); }\n")
					f.write("\n")
	f.write("#endif // " + header_name + "\n")
