#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <exception>
#include <chrono>
#include <sstream>

#include <cmath>
#include <cstdlib>

class Expression {
    enum BINARY_OPERATOR {
        OPERATOR_SUM,
        OPERATOR_DIF,
        OPERATOR_PROD,
        OPERATOR_DIV,
        OPERATOR_POW
    };
    // общие данные для всех экземпляров классов
    // необходимы для универсальности класса Expression и 
    // простой расширяемости: добавил функцию и имя в список функций и забыл об этом
    static std::vector<double>      vars;               // Переменные (необходимость данного списка пока под вопросом)       
    static std::vector<std::string> varSymbols;         // Имена переменных
    static std::vector<std::string>        funcUName;   // Имена известных функций 1 переменной
    static std::vector<double (*)(double)> funcU;       // Сами функции 1 переменной
    static std::vector<std::string>                funcBName; // Предполагалось для хранения знаков бинарных операций
    static std::vector<double (*)(double, double)> funcB;     // Сами бинарные операции
    static std::vector<int>                        funcBLevel;
    static int maxLevel;
    // Необъявленные явно функции
    static double  neg(double a);
    static double  sum(double a, double b);
    static double  dif(double a, double b);
    static double prod(double a, double b);
    static double  div(double a, double b);
    // Вспомогательные функции
    static int check_brackets(std::string a);
    static bool is_number(std::string line);
public:
    // Это абстрактный класс, который ничего сам не хранит, кроме общих данных
    virtual ~Expression();
    // Основные функции - 
    virtual double Calc() = 0;
    virtual Expression Calc(std::vector<std::tuple<std::string, Expression*>>) = 0;
    virtual std::string Print(int level = 0) = 0;
    // Генерация деревьев выражений
    static std::string PrepareR(std::string expr);
    static std::string Prepare(std::string);
    static Expression* ParseExpressionR(std::string expr);
    static Expression* ParseExpression(std::string expr);
    // Функции для работы с общими данными
    // с переменными
    static void SetVariables(std::vector<double> v, std::vector<std::string> s);
    static void PrintVariables();
    static double GetVariable(size_t index);
    static std::string GetVariableName(size_t index);
    static size_t AddVariable(std::string varSymbol);
    // с функциями 1 переменной
    static size_t GetNFuncU();
    static double GetFuncU(size_t index, double x);
    static std::string GetFuncUName(size_t index);
    static size_t AddFunctionU(std::string name, double (*func)(double));
    // с функциями 2 переменных
    inline static int GetMaxLevel();
    static int GetFuncBLevel(size_t index);
    static size_t GetNFuncB();
    static double GetFuncB(size_t index, double x, double y);
    static std::string GetFuncBName(size_t index);
    static size_t AddFunctionB(std::string name, double (*func)(double, double), int level);
    // Инициализация базовых переменных, требуется выполнять один раз перед использованием всех классов
    static void Init();
};

class TermConstDouble : public Expression {
public:
    double arg;

    TermConstDouble();
    TermConstDouble(double arg);
    ~TermConstDouble();

    double Calc();
    std::string Print(int level = 0);

};

class TermVariable : public Expression {
public:
    int index;

    TermVariable();
    TermVariable(int index);
    TermVariable(std::string var);
    ~TermVariable();

    double Calc();
    std::string Print(int level = 0);
};

class TermBinary: public Expression {
public:
    Expression *left;
    Expression *right;
    size_t operation; // + - * /

    TermBinary(Expression *left, Expression *right, size_t operation);
    ~TermBinary();

    double Calc();
    std::string Print(int level = 0);
};

class TermUnary : public Expression {
public:
    Expression *next;
    size_t operation;

    TermUnary(Expression *next, size_t operation);
    ~TermUnary();

    double Calc();
    std::string Print(int level = 0);
};

#ifdef HAS_TEST

double test(double x, double y, double z);
void test_expression();

#endif // HAS_TEST