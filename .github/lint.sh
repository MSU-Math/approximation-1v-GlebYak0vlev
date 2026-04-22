#!/usr/bin/env bash
ERRORS=0
function countError {
    echo "Проверка завершилась с ошибкой."
    (( ERRORS++ ))
}
trap countError ERR 

make clean || true

echo "# cppcheck"
cppcheck --enable=warning,performance --error-exitcode=2 . && echo "OK" || false

echo "# clang-format (форматирование кода)"
find . -maxdepth 1 -name "*.c*" | xargs -r $(dirname $0)/clang-format.sh --dry-run -Werror && echo "OK" || false
find . -maxdepth 1 -name "*.h" | xargs -r $(dirname $0)/clang-format.sh --dry-run -Werror && echo "OK" || false

echo "# Предупреждения GCC"
$(dirname $0)/build.sh && echo "OK" || false

exit $ERRORS
