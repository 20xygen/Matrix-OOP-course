# Matrix

## Техническое задание

Написать шаблонный класс Residue<size_t N> – кольцо вычетов по модулю N (которое, как известно, в случае простого N является полем). Должны поддерживаться арифметические операции, кроме деления, а в случае простого N и деление тоже. Элемент кольца Residue<N> должно быть можно сконструировать от int путем взятия остатка от деления на N (в математическом смысле). Также должно быть можно привести Residue<N> к int (но преобразования должны быть явными в обе стороны). Попытка вызвать операцию деления при составном N должна приводить к ошибке компиляции. Проверка числа N на простоту должна выполняться в compile-time, причем время работы должно быть O(sqrt(N)).

После этого, используя ранее написанный класс рациональных чисел, написать класс Matrix с тремя шаблонными параметрами: size_t M, size_t N, typename Field=Rational. (По умолчанию берётся поле рациональных чисел, но можно создать матрицу и над другим полем.)

Конструктор матрицы по умолчанию должен возвращать матрицу из одних нулей. Можно добавить статический метод unityMatrix(), возвращающий единичную матрицу (но если M != N, это не должно компилироваться).

Матрицы должны поддерживать следующие операции:

- Сложение, вычитание, операторы +=, -=. Сложение и вычитание матриц несоответствующих размеров не должно компилироваться.
- Умножение на число типа Field.
- Умножение двух матриц. Попытка перемножить матрицы несоответствующих размеров должна приводить к ошибке компиляции.
- Метод det(), возвращающий определитель матрицы за O(N3). Взятие определителя матрицы, у которой M != N, не должно компилироваться. Определитель матрицы над кольцом, которое не является полем, можно не поддерживать.
- Метод transposed(), возвращающий транспонированную матрицу.
- Метод rank() - вычислить ранг матрицы. Ранг матрицы над кольцом, которое не является полем, можно не поддерживать. Вычисление ранга должно работать за O(N3).
- Методы inverted() и invert() - вернуть обратную матрицу и обратить данную матрицу. Попытка обратить матрицу, у которой M != N, должна приводить к ошибке компиляции.
- Метод trace() - вычислить след матрицы. Вычисление следа от неквадратной матрицы не должно компилироваться.
- Методы getRow(unsigned) и getColumn(unsigned), возвращающие std::array из соответствующих значений.
- К матрице должен быть дважды применим оператор [], причём это должно работать как для неконстантных, так и для константных матриц. В первом случае содержимое матрицы должно быть можно таким способом поменять.

Другие способы изменения содержимого матрицы, кроме описанных выше, должны отсутствовать. Однако не запрещается реализовать дополнительные методы для выполнения каких-либо иных алгебраических операций или для удобства работы, если по названию и сигнатуре этих методов будет без комментариев понятно их действие. Квадратные матрицы размера N должно быть можно объявлять всего с одним обязательным шаблонным параметром: SquareMatrix<N>.

В вашем файле должна отсуствовать функция main, а сам файл должен называться matrix.h. Ваш код будет вставлен посредством #include в программу, содержащую тесты. В качестве компилятора необходимо выбирать Make.
