from setuptools import setup, Extension

module_test = Extension('impf.clib_test',
			sources=['src/pytest.c'])
module_optm = Extension('impf.clib_optm',
			sources=['src/pyoptm.c',
	    			'src/utils.c',
				'src/linalg/daxpy.c',
				'src/linalg/dscal.c',
				'src/lp/simplex_std.c',
				'src/lp/simplex_gen.c'],
			include_dirs=['include'])

setup(
	name='impf',
	version='1.0',
	author='Zhuang Linsheng',
	author_email='zhuanglinsheng@outlook.com',
	description='Python interface for impf C library',
	license='GPL3',
	url='https://github.com/zhuanglinsheng/impf',
	packages=['impf'],
	ext_modules=[
		module_test,
		module_optm,
		],
	python_requires=">= 3.8",
	classifiers=[
		'Programming Language :: Python :: 3',
		'License :: OSI Approved :: GPL3 License',
		'Operating System :: OS Independent',
	],
)
