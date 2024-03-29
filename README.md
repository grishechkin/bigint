# bigint

## Реализация

Класс называется `big_integer`, а код решения находится в файлах `big_integer.h` и `big_integer.cpp`.

Реализованны:
- Конструктор по умолчанию, инициализирующий число нулём.
- Конструктор копирования.
- Конструкторы от числовых типов.
- Explicit конструктор от `std::string`.
- Оператор присваивания.
- Операторы сравнения.
- Арифметические операции: сложение, вычитание, умножение, деление, унарный минус и плюс.
- Инкреметы и декременты.
- Битовые операции: и, или, исключающее или, не (аналогично битовым операциям для `int`)
- Битовые сдвиги.
- Внешняя функция `std::string to_string(big_integer const&)`.

Реализация удовлетворяет следующим требованиям:
- Умножение и деление должны работать не хуже, чем за O(nm).
- Остальные операции должны производиться с максимально возможной асимптотической эффективностью.
- Помимо асимптотики, стоит уделить внимание оптимизации количества аллокаций и общего времени работы.
- Разряды числа должны представляться как минимум 32-битными числами, все биты в их представлении должны использоваться.
- `big_integer` должен уметь создаваться от числовых типов сохраняя значение. Если переменная числового типа имела значение `x`, значит `big_integer` после создания должен иметь значение `x`.
- Время прохождения тестов в CI не должно превышать 1 секунду в Release-сборке, при этом ваше решение должно укладываться в лимит при каждом перезапуске.

