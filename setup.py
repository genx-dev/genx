from setuptools import setup, Extension
import genx.version


def build_ext():
    try:
        from Cython.Build import cythonize
    except ImportError:
        print('Cython is not available, cython extension not compiled.')
        extentions = list()
    else:
        #kwargs = {'extra_link_args': ['-Wl,--strip-all'], 'extra_compile_args': ['-Wno-cpp']}
        kwargs = {'extra_link_args': ['-Wl'], 'extra_compile_args': ['-Wno-cpp']}
        offspec = Extension("genx.models.lib.offspec_ext", ["genx/models/lib/offspec_ext.pyx"], **kwargs)
        paratt = Extension("genx.models.lib.paratt_ext", ["genx/models/lib/paratt_ext.pyx"], **kwargs)
        sxrd = Extension("genx.models.lib.sxrd_ext", ["genx/models/lib/sxrd_ext.pyx"], **kwargs)
        return cythonize([offspec, paratt, sxrd], nthreads=4, compiler_directives={'language_level': 3, 'cdivision': True, 'boundscheck': False, 'wraparound': False})


_data_files = [('/usr/share/doc/genx/', ['changelog.txt', 'README.txt']), ('/etc/genx/', ['genx/genx.conf']), ('/etc/genx/profiles/', ['genx/profiles/Default.conf', 'genx/profiles/Reflectivity.conf', 'genx/profiles/SuperAdam.conf']), ('/usr/share/licenses/genx/', ['LICENSE.txt']), ('/usr/share/applications/', ['resources/genx.desktop']), ('/usr/share/pixmaps/', ['resources/genx.xpm'])]

for i in ('16', '22', '32', '48', '64', '128'):
   _data_files.append(('/usr/share/icons/hicolor/' + i + 'x' + i + '/', ['resources/genx_' + i + 'x' + i + '.png']))



setup(
    name="GenX",
    version=genx.version.version,
    description="X-ray and Neutron reflectivity fitting software.",
    scripts=['genx/genx'],
    packages=['genx','genx.plugins', 'genx.plugins.add_ons',  'genx.plugins.add_ons.help_modules',
        'genx.plugins.data_loaders', 'genx.plugins.data_loaders.help_modules',
        'genx.models', 'genx.models.lib', 'genx.lib'],
    package_data={
                  'genx': ['examples/*.*'], 
                  'genx.models': ['databases/*.*', 'databases/f1f2_cxro/*.*', 'databases/f1f2_nist/*.*'], 
                  },
    data_files = _data_files,
    ext_modules=build_ext(),
    author='Matts Bjorck, ported to Python 3 by Hao Zhang',
    license='GPL v3',
    author_email = "matts.bjorck@gmail.com",
    url = "http://genx.sourceforge.net/",
    install_requires=['numpy', 'scipy', 'matplotlib', 'appdirs', 'wxpython', 'appdirs', 'h5py'],
)
