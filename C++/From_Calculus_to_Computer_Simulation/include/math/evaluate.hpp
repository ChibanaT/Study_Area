#ifndef EVALUATE_EXPR_HPP
#define EVALUATE_EXPR_HPP

#include <string>
#include <stack>
#include <queue>
#include <vector>
#include <unordered_map>
#include <functional>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace expr {

// Helper: check if a character is operator
inline bool isOperator(char c) {
    return c == '+' || c == '-' || c == '*' || c == '/' || c == '^';
}

// Operator precedence
inline int precedence(char op) {
    switch (op) {
        case '+': case '-': return 1;
        case '*': case '/': return 2;
        case '^': return 3;
    }
    return 0;
}

// Operator associativity (true = right-associative)
inline bool isRightAssociative(char op) {
    return (op == '^');
}

// --- Convert infix expression to RPN using Shunting Yard ---
inline std::queue<std::string> toRPN(const std::string& expr, double x) {
    std::queue<std::string> output;
    std::stack<std::string> ops;

    std::unordered_map<std::string, int> precedenceMap = {
        {"+", 1}, {"-", 1}, {"*", 2}, {"/", 2}, {"^", 3}
    };

    std::unordered_map<std::string, bool> rightAssoc = {
        {"^", true}
    };

    for (size_t i = 0; i < expr.size();) {
        if (std::isspace(expr[i])) {
            ++i;
            continue;
        }

        // --- Number or variable ---
        if (std::isdigit(expr[i]) || expr[i] == '.') {
            std::string num;
            while (i < expr.size() && (std::isdigit(expr[i]) || expr[i] == '.'))
                num += expr[i++];
            output.push(num);
            continue;
        }

        // --- Handle std::function_name or function_name ---
        if (std::isalpha(expr[i]) || 
           (i + 4 < expr.size() && expr.substr(i, 5) == "std::")) {
            std::string name;
            
            // Skip "std::" if present
            if (i + 4 < expr.size() && expr.substr(i, 5) == "std::") {
                i += 5; // Skip "std::"
            }
            
            // Read the function/variable name
            while (i < expr.size() && (std::isalpha(expr[i]) || expr[i] == '_'))
                name += expr[i++];

            if (name == "x") {
                output.push(std::to_string(x));
            } else {
                // function (sin, cos, log, pow, etc)
                ops.push(name);
            }
            continue;
        }

        // --- Handle comma (for multi-argument functions) ---
        if (expr[i] == ',') {
            while (!ops.empty() && ops.top() != "(") {
                output.push(ops.top());
                ops.pop();
            }
            ++i;
            continue;
        }

        // --- Operator ---
        if (isOperator(expr[i])) {
            std::string op(1, expr[i]);
            while (!ops.empty() && 
                   precedenceMap.count(ops.top()) &&
                   ((rightAssoc.count(op) == 0 && precedenceMap[op] <= precedenceMap[ops.top()]) ||
                    (rightAssoc.count(op) && precedenceMap[op] < precedenceMap[ops.top()]))) {
                output.push(ops.top());
                ops.pop();
            }
            ops.push(op);
            ++i;
            continue;
        }

        // --- Parentheses ---
        if (expr[i] == '(') {
            ops.push("(");
            ++i;
            continue;
        }

        if (expr[i] == ')') {
            while (!ops.empty() && ops.top() != "(") {
                output.push(ops.top());
                ops.pop();
            }
            if (!ops.empty() && ops.top() == "(") {
                ops.pop();
            }
            // if function before '('
            if (!ops.empty() && std::isalpha(ops.top()[0])) {
                output.push(ops.top());
                ops.pop();
            }
            ++i;
            continue;
        }

        throw std::runtime_error("Unexpected character in expression");
    }

    while (!ops.empty()) {
        output.push(ops.top());
        ops.pop();
    }

    return output;
}

// --- Evaluate RPN ---
inline double evalRPN(std::queue<std::string> rpn) {
    std::stack<double> st;

    std::unordered_map<std::string, std::function<double(double)>> funcs1 = {
        {"sin", [](double a){ return std::sin(a); }},
        {"cos", [](double a){ return std::cos(a); }},
        {"tan", [](double a){ return std::tan(a); }},
        {"exp", [](double a){ return std::exp(a); }},
        {"log", [](double a){ return std::log(a); }},
        {"sqrt", [](double a){ return std::sqrt(a); }},
        {"abs", [](double a){ return std::fabs(a); }}
    };

    std::unordered_map<std::string, std::function<double(double,double)>> funcs2 = {
        {"+", [](double a, double b){ return a + b; }},
        {"-", [](double a, double b){ return a - b; }},
        {"*", [](double a, double b){ return a * b; }},
        {"/", [](double a, double b){ return a / b; }},
        {"^", [](double a, double b){ return std::pow(a, b); }},
        {"pow", [](double a, double b){ return std::pow(a, b); }}  // Added pow function
    };

    while (!rpn.empty()) {
        std::string token = rpn.front();
        rpn.pop();

        if (funcs2.count(token)) {
            if (st.size() < 2) {
                throw std::runtime_error("Not enough operands for binary operation");
            }
            double b = st.top(); st.pop();
            double a = st.top(); st.pop();
            st.push(funcs2[token](a, b));
        }
        else if (funcs1.count(token)) {
            if (st.empty()) {
                throw std::runtime_error("Not enough operands for unary operation");
            }
            double a = st.top(); st.pop();
            st.push(funcs1[token](a));
        }
        else {
            // number
            try {
                st.push(std::stod(token));
            } catch (const std::exception&) {
                throw std::runtime_error("Invalid number: " + token);
            }
        }
    }

    if (st.size() != 1) {
        throw std::runtime_error("Invalid expression: stack size is " + std::to_string(st.size()));
    }

    return st.top();
}

// --- Main evaluate function ---
inline double evaluate(const std::string& expr, double x) {
    auto rpn = toRPN(expr, x);
    return evalRPN(rpn);
}

} // namespace expr

#endif // EVALUATE_EXPR_HPP