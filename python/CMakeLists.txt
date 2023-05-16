
execute_process(
	COMMAND "${PYTHON_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_paths()['purelib'])"
	OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
	OUTPUT_STRIP_TRAILING_WHITESPACE
	)

# install(TARGETS paint
#	COMPONENT python
#	LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}"
#	)

