#include <iostream>
#include "num_methodsp2.h"
#include "nums_tests.h"


// Просто расскомментируйте метод, который хотите протестировать
// Тест можно менять, передав параметры, отличные от аргументов по умолчанию
// Сами численные методы лежат в "num_methodsp2.h"
// Внимательно читайте описание к каждой функции! Какие-то нужны только для тестирования

int main() {
    setlocale(LC_ALL, "Russian");
    //nums::test_ERK2();
    //nums::test_ROS_1_a05();
    nums::test_ROS_1_a1();
    //nums::test_ERK_2_prec();
    //nums::TestShootingMethod();
    //nums::TestGridMethod();
    return 0;
}
