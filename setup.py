from distutils.core import setup
import py2exe
##TMP
import matplotlib

setup(
    console=['gui_app.py'],
    options={
             'py2exe': {
                        'packages' : ['matplotlib', 'pytz'],
                        'excludes': ['_gtkagg', '_tkagg', '_agg2', '_cairo', '_cocoaagg',
                                     '_fltkagg', '_gtk', '_gtkcairo', '_tkagg'],
                       }
            },
    data_files=matplotlib.get_py2exe_datafiles()
)
