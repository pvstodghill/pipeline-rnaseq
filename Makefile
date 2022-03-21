default: doit09.out

clean :
	rm -rf data doit*.out

doit%.out : doit-preamble.bash config.bash

doit00.out: doit00.bash
	run ./doit00.bash

doit01.out: doit01.bash doit00.out
	run ./doit01.bash


doit02.out: doit02.bash doit01.out
	run ./doit02.bash


doit03.out: doit03.bash doit02.out
	run ./doit03.bash


doit04.out: doit04.bash doit03.out
	run ./doit04.bash


doit05.out: doit05.bash doit04.out
	run ./doit05.bash


doit06.out: doit06.bash doit05.out
	run ./doit06.bash


doit07.out: doit07.bash doit06.out
	run ./doit07.bash


doit08.out: doit08.bash doit07.out
	run ./doit08.bash


doit09.out: doit09.bash doit08.out
	run ./doit09.bash

