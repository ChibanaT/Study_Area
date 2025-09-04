#ifndef EXPR_TO_CMATH_HPP
#define EXPR_TO_CMATH_HPP

#include <string>
#include <regex>
#include <unordered_map>

namespace expr {

// Function to remove all whitespace characters
inline std::string removeWhitespace(const std::string& str) {
    std::string result;
    result.reserve(str.length()); // Reserve space for efficiency
    
    for (size_t i = 0; i < str.length(); ++i) {
        char c = str[i];
        // Skip all whitespace characters
        if (c != ' ' && c != '\t' && c != '\n' && c != '\r' && c != '\v' && c != '\f') {
            result += c;
        }
    }
    
    return result;
}

inline std::string to_cmath(const std::string& expr) {
    // 0. Remove ALL whitespace first
    std::string res = removeWhitespace(expr);

    // 1. Replace mathematical functions
    std::unordered_map<std::string, std::string> funcs = {
        {"sin", "std::sin"},
        {"cos", "std::cos"},
        {"tan", "std::tan"},
        {"exp", "std::exp"},
        {"log", "std::log"}, // natural log
        {"sqrt", "std::sqrt"},
        {"abs", "std::abs"},
        {"floor", "std::floor"},
        {"ceil", "std::ceil"}
    };
    for (auto& [k, v] : funcs) {
        res = std::regex_replace(res, std::regex("\\b" + k + "\\b"), v);
    }

    // 2. Replace pi constant
    res = std::regex_replace(res, std::regex("\\bpi\\b"), "3.14159265358979323846");

    // 3. Handle e^something first (without spaces now)
    // e^4 -> std::exp(4)
    std::regex e_pow_regex(R"(\be\^(\([^()]+\)|[a-zA-Z0-9\+\-\*/]+))");
    while (std::regex_search(res, e_pow_regex)) {
        res = std::regex_replace(res, e_pow_regex, "std::exp($1)");
    }

    // 4. Replace standalone e
    res = std::regex_replace(res, std::regex("\\be\\b"), "std::exp(1.0)");

    // 5. Handle powers: a^b -> std::pow(a,b) (without spaces)
    std::regex pow_regex(R"((\([^\)]+\)|[a-zA-Z]+|[0-9]+)\^(\([^\)]+\)|-?[a-zA-Z0-9]+))");
    while (std::regex_search(res, pow_regex)) {
        res = std::regex_replace(res, pow_regex, "std::pow($1,$2)");
    }

    // 6. Insert multiplication for implicit multiplication: 2x -> 2*x, 2(x+1) -> 2*(x+1)
    std::regex mult_regex(R"((\d)([a-zA-Z\(]))");
    res = std::regex_replace(res, mult_regex, "$1*$2");

    // 7. Insert multiplication between parentheses: (2+1)(3+2) -> (2+1)*(3+2)
    std::regex paren_mult_regex(R"((\))(\())");
    res = std::regex_replace(res, paren_mult_regex, "$1*$2");

    return res;
}

} // namespace expr

#endif // EXPR_TO_CMATH_HPP