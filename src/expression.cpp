#include "expression.h"

// Общие данные
std::vector<double>                     Expression::vars       ;
std::vector<std::string>                Expression::varSymbols ;
std::vector<std::string>                Expression::funcUName  ;
std::vector<double (*)(double)>         Expression::funcU      ;
std::vector<std::string>                Expression::funcBName  ;
std::vector<double (*)(double, double)> Expression::funcB      ;
std::vector<int>                        Expression::funcBLevel ;
int Expression::maxLevel;

Expression::~Expression() {};
// Отсутствующие функции-операторы
double Expression::neg(double a) { return -a; };
double Expression::sum(double a, double b)  { return a + b; };
double Expression::dif(double a, double b)  { return a - b; };
double Expression::prod(double a, double b) { return a * b; };
double Expression::div(double a, double b)  { return a / b; };

// Работа с переменными
void Expression::SetVariables(std::vector<double> v, std::vector<std::string> s) {
    if (v.size() == s.size()) {
        int has_var = 0;
        for (size_t j = 0; j< v.size(); j++) {
            for (size_t i = 0; i<varSymbols.size(); i++) {
                if (varSymbols[i] == s[j]) {
                    vars[i] = v[j];
                    has_var = 1;
                    break;
                }
            }
            if (has_var == 0) {
                varSymbols.push_back(s[j]);
                vars.push_back(v[j]);
            }
        }
    }
}

void Expression::PrintVariables() {
    for (size_t i = 0; i<varSymbols.size(); i++) {
        std::cout << varSymbols[i] << " = " << vars[i] << "; ";
    }
    std::cout << std::endl;
}

double Expression::GetVariable(size_t index) {
    return vars[index];
}

std::string Expression::GetVariableName(size_t index) {
    return varSymbols[index];
}

size_t Expression::AddVariable(std::string varSymbol) {
    for (size_t i = 0; i<varSymbols.size(); i++) {
        if (varSymbols[i] == varSymbol)
            return (int) i;
    }
    varSymbols.push_back(varSymbol);
    vars.push_back(0);
    return varSymbols.size() - 1;
}

// с функциями 1 переменной
size_t Expression::GetNFuncU() {
    return funcU.size();
}

double Expression::GetFuncU(size_t index, double x) {
    return funcU[index](x);
}

std::string Expression::GetFuncUName(size_t index) {
    return funcUName[index];
}

size_t Expression::AddFunctionU(std::string name, double (*func)(double)) {
    for (size_t i = 0; i<funcUName.size(); i++) {
        if (funcUName[i] == name)
            return (int) i;
    }
    funcUName.push_back(name);
    funcU.push_back(func);
    return varSymbols.size() - 1;
}

// с функциями 2 переменных
inline int Expression::GetMaxLevel() {
    return maxLevel;
}

int Expression::GetFuncBLevel(size_t index) {
    return funcBLevel[index];
}

size_t Expression::GetNFuncB() {
    return funcB.size();
}

double Expression::GetFuncB(size_t index, double x, double y) {
    return funcB[index](x, y);
}

std::string Expression::GetFuncBName(size_t index) {
    return funcBName[index];
}

size_t Expression::AddFunctionB(std::string name, double (*func)(double, double), int level) {
    for (size_t i = 0; i<funcBName.size(); i++) {
        if (funcBName[i] == name)
            return (int) i;
    }
    funcBName.push_back(name);
    funcB.push_back(func);
    funcBLevel.push_back(level);
    return varSymbols.size() - 1;
}

// Инициализация базовых переменных, требуется выполнять один раз перед использованием всех классов
void Expression::Init() {
    AddFunctionU("-"   , neg);
    AddFunctionU("exp" , exp);
    AddFunctionU("log" , log);
    AddFunctionU("sin" , sin);
    AddFunctionU("cos" , cos);
    AddFunctionU("tan" , tan);
    AddFunctionU("atan", atan);
    AddFunctionU("sinh", sinh);
    AddFunctionU("cosh", cosh);
    AddFunctionU("tanh", tanh);
    AddFunctionU("erf" , erf);

    AddFunctionB("+", sum , 0);
    AddFunctionB("-", dif , 0);
    AddFunctionB("*", prod, 1);
    AddFunctionB("/", div , 1);
    AddFunctionB("^", pow , 2);
    maxLevel = 2;
    //AddFunctionB("atan2", atan2);
}

// Константа
TermConstDouble::TermConstDouble() : arg(0) {};
TermConstDouble::TermConstDouble(double arg) : arg(arg) {};
TermConstDouble::~TermConstDouble() {};

double TermConstDouble::Calc() { 
    return arg; 
};

std::string TermConstDouble::Print(int level) {
    std::ostringstream ss;
    ss << arg;
    return ss.str(); 
};

// Переменная
TermVariable::TermVariable() : index(-1) {};
TermVariable::TermVariable(int index) : index(index) {};
TermVariable::TermVariable(std::string var) {
    AddVariable(var);
};
TermVariable::~TermVariable() {};

double TermVariable::Calc() { 
    return (index != -1) ? GetVariable(index) : 0; 
};

std::string TermVariable::Print(int level) { 
    return (index != -1) ? GetVariableName(index) : ""; 
};

// Бинарная операция
TermBinary::TermBinary(Expression *left, Expression *right, size_t operation) : 
    left(left), right(right), operation(operation) {
};

TermBinary::~TermBinary() {
    delete left;
    delete right;
    left  = nullptr;
    right = nullptr;
};

double TermBinary::Calc() { 
    return GetFuncB(operation, left->Calc(), right->Calc());
};

std::string TermBinary::Print(int level) {
    return GetFuncBName(operation) + "(" + left->Print() + "," + right->Print() + ")";
}

// функция 1 переменной
TermUnary::TermUnary(Expression *next, size_t operation) : next(next), operation(operation) {};
TermUnary::~TermUnary() {
    delete next;
    next = nullptr;
};

double TermUnary::Calc() {
    try { 
        return GetFuncU(operation, next->Calc());
    } catch (std::exception &) {
        std::cout << "Unknown operation" << std::endl;
        throw std::exception();
    }
};

std::string TermUnary::Print(int level) {
    return GetFuncUName(operation) + "(" + next->Print() + ")";
}

// Вспомогательные функции класса
int Expression::check_brackets(std::string a) {
    int bracket = 0;
    for (size_t i = 0; i<a.size(); i++) {
        if (a[i] == '(') 
            bracket++;
        if (a[i] == ')')
            bracket--;
        if (bracket < 0) {
            std::cout << "Bad expression: check brackets in position: " << i << std::endl;
            return i; 
        }
    }
    return bracket;
}

bool Expression::is_number(std::string line)
{
    char* p;
    strtod(line.c_str(), &p);
    return *p == 0;
}

std::string Expression::PrepareR(std::string expr) {
    std::string pre, term1, term2;
    int level = 3, maxLevel = 3;
    size_t size = expr.size();
    int bracket = 0, brackets = 0;
    for (size_t i = 0; i<size; i++) {
        if (expr[i] == '(') {
            if (bracket == 0) {
                brackets++;
            }
            bracket++;
        }
        if (expr[i] == ')') {
            bracket--;
            if ((i == size - 1) && (brackets == 1) && (expr[0] == '(')){
                expr.erase(0, 1);
                expr.erase(size - 2, 1);
            } 
        }
    }
    
    size = expr.size();
    for (size_t i = 0; i<size; i++) {
        if (expr[i] == '(') {
            if (bracket == 0) {
                brackets++;
                pre.assign(expr, 0, i);
                term1.assign(expr, i, size - i);
            }
            bracket++;
        }
        if (expr[i] == ')') {
            bracket--;
        }
    }

    // std::cout << "expr: " << expr << std::endl;
    
    const char *op = "+-*/^";
    int opl[5] = {0, 0, 1, 1, 2};
    int opi;

    for (size_t i = 0; i<size; i++) {
        if (expr[i] == '(') {
            bracket++;
        }
        if (expr[i] == ')') {
            bracket--;
        }
        for (size_t j = 0; j<5; j++) {
            if (expr[i] == op[j]) {
                if ((bracket == 0) && (level > opl[j])) {
                    term1.assign(expr, 0, i);
                    term2.assign(expr, i + 1, size - i);
                    level = opl[j];
                    opi = j;
                }
            }
        }
    }
    if (level < maxLevel) {
        std::string str;
        str += op[opi];
        return str + "(" + PrepareR(term1) + "," + PrepareR(term2) + ")";
    }
    // std::cout << "pre: " << pre << std::endl;
    // std::cout << "term1: " << term1 << std::endl;
    // std::cout << "term2: " << term2 << std::endl;
    // getchar();
    if (brackets > 0) {
        return pre + "(" + PrepareR(term1) + ")";
    } else {
        return expr;
    }
}

std::string Expression::Prepare(std::string expr) {
    std::size_t k = expr.find(' ');
    while ( k != std::string::npos) {
        expr.erase(k, 1);
        k = expr.find(' ');
    };

    return PrepareR(expr);        
}

// Рекурсивный парсер выражения без пробелов
Expression* Expression::ParseExpressionR(std::string expr) {
    int bracket = 0, brackets = 0;
    std::string pre, term1, term2;
    size_t size = expr.size();
    
    for (size_t i = 0; i<size; i++) {
        if (expr[i] == '(') {
            if (bracket == 0) {
                brackets++;
            }
            bracket++;
        }
        if (expr[i] == ')') {
            bracket--;
            if ((i == size - 1) && (brackets == 1) && (expr[0] == '(')){
                expr.erase(0, 1);
                expr.erase(size - 2, 1);
            } 
        }
    }

    for (size_t i = 0; i<size; i++) {
        if (expr[i] == '(') {
            if (bracket == 0) {
                pre.assign(expr, 0, i);
                term1.assign(expr, i, size - i);
            }
            bracket++;
        }
        if (expr[i] == ')') {
            bracket--;
        }
    }
    // std::cout << "after brackets deletion: " << expr << std::endl;

    int maxLevel = Expression::GetMaxLevel() + 1; 
    int level = maxLevel;
    size_t op;

    size = expr.size();
    for (size_t i = 0; i<size; i++) {
        if (expr[i] == '(') {
            bracket++;
        }
        if (expr[i] == ')') {
            bracket--;
        }
        for (size_t j = 0; j<Expression::GetNFuncB(); j++) {
            if (expr.compare(i, Expression::GetFuncBName(j).size(), Expression::GetFuncBName(j)) == 0) {
                if ((bracket == 0) && (level > Expression::GetFuncBLevel(j))) {
                    term1.assign(expr, 0, i);
                    term2.assign(expr, i + 1, size - i);
                    level = Expression::GetFuncBLevel(j);
                    op = j;
                    // std::cout << "pre: " << pre << std::endl;
                    // std::cout << "term1: " << term1 << std::endl;
                    // std::cout << "term2: " << term2 << std::endl;
                }
            }
        }
    }
    // std::cout << "pre: " << pre << std::endl;
    // std::cout << "term1: " << term1 << std::endl;
    // std::cout << "term2: " << term2 << std::endl;
    // std::cout << "next block: " << std::endl;
    if (level < maxLevel) {
        return new TermBinary(ParseExpressionR(term1), ParseExpressionR(term2), op);
    }

    for (size_t j = 0; j<Expression::GetNFuncU(); j++) {
        if (pre.compare(Expression::GetFuncUName(j)) == 0) {
            return new TermUnary(ParseExpressionR(term1),  j);
        }
    }

    std::cout << "const: " << expr << " " << atof(expr.c_str()) << std::endl;
    if (is_number(expr)) {
        // std::cout << "const: " << expr << " " << atof(expr.c_str()) << std::endl;
        return new TermConstDouble(atof(expr.c_str()));
    } else {
        // std::cout << "variable: " << expr << std::endl;
        size_t index = Expression::AddVariable(expr);
        return new TermVariable(index);
    }
}

// Точка вызова рекурсивного парсера для выражения с пробелами
Expression* Expression::ParseExpression(std::string expr) {
    if (expr.empty()) 
        return new TermConstDouble();
    std::cout << expr << std::endl;
    int k = 0;
    std::size_t i = expr.find(' ');
    while ( i != std::string::npos) {
        // std::cout << k << ": " << i << " " << expr << std::endl;
        expr.erase(i, 1);
        k++;
        i = expr.find(' ');
    };

    // std::cout << expr << std::endl;
    if (check_brackets(expr))
        return new TermConstDouble();

    return ParseExpressionR(expr);
}

#ifdef HAS_TEST

double test(double x, double y, double z) {
    return sin(pow(x*y + y*z + exp(x*(z - 1)), 1.5));
}

void test_expression() {
    Expression::Init();
    Expression *expr;
    expr = Expression::ParseExpression("exp(-x^2)");
    std::cout << "result: " << expr->Print() << std::endl;
    expr = Expression::ParseExpression("(x + a)*exp(-x^2)");
    std::cout << "result: " << expr->Print() << std::endl;
    expr = Expression::ParseExpression("((x + a)*(y - b))");
    std::cout << "result: " << expr->Print() << std::endl;
    expr = Expression::ParseExpression("sin(exp(-x^2))");
    std::cout << "result: " << expr->Print() << std::endl;
    expr = Expression::ParseExpression("cos(exp(-x^2)*pi*cos(4*y))");
    std::cout << "result: " << expr->Print() << std::endl;
    expr = Expression::ParseExpression("sin((x*y + y*z + exp(x*(z - 1)))^1.5)");
    std::cout << "result: " << expr->Print() << std::endl;
    std::cout << "parser: " << Expression::Prepare("sin((x*y + y*z + exp(x*(z - 1)))^1.5)") << std::endl;
    expr->SetVariables({1, 2, 3}, {"x", "y", "z"});
    expr->PrintVariables();

    using namespace std::chrono;

    double result;
    steady_clock::time_point t1, t2;
    duration<double> time_span, time_span0;
    int Nt = 1000000;

    t1 = steady_clock::now();
    for (int i = 0; i<Nt; i++) {
        result = 1;
    }
    t2 = steady_clock::now();
    time_span0 = duration_cast<duration<double>>(t2 - t1);

    t1 = steady_clock::now();
    for (int i = 0; i<Nt; i++) {
        result = 1;
        result = result*test(1, 2, 3);
    }
    t2 = steady_clock::now();
    time_span = duration_cast<duration<double>>(t2 - t1);
    std::cout << "CPP result: " << result << " time: " << (time_span.count() - time_span0.count())/Nt << std::endl;

    t1 = steady_clock::now();
    for (int i = 0; i<Nt; i++) {
        result = 1;
        result = result*expr->Calc();
    }
    t2 = steady_clock::now();
    time_span = duration_cast<duration<double>>(t2 - t1);
    std::cout << "my  result: " << result << " time: " << (time_span.count() - time_span0.count())/Nt << std::endl;
    delete expr;
}

#endif // HAS_TEST

