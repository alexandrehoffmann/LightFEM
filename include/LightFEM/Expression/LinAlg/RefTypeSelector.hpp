/*
 * RefTypeSelector.hpp
 * 
 * Copyright 2022 Alexandre Hoffmann <alexandre.hoffmann.etu@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#ifndef REF_TYPE_SELECTOR_HPP
#define REF_TYPE_SELECTOR_HPP

#include <type_traits>
#include <utility>

template < typename Expr >
class RefTypeSelector
{
private:
	using Expr1 = typename std::decay<Expr>::type;
public:
	using Type = typename std::conditional< std::is_lvalue_reference<Expr>::value, const Expr1&, Expr1>::type;
};

#endif // REF_TYPE_SELECTOR_HPP
