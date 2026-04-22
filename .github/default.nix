with import <nixpkgs> {};

mkShell {
    buildInputs = [ gcc gdb 
        cppcheck clangAnalyzer clang-tools
        qt5.qtbase qt5.qmake qt5.qtwayland
        python3Packages.compiledb
    ];
}
