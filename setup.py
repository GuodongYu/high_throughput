from setuptools import setup,find_packages

NAME='high_throughput'
PACKAGES = [NAME] + ["%s.%s" % (NAME, i) for i in find_packages(NAME)]

setup (
        name = NAME,
        version = '0.0.1',
        packages = PACKAGES,
        author = 'Guodong Yu',
        description = 'high_throughput',
        zip_safe = False,
        long_description = 'Often used hight throughput tools for Guodong',
        author_email = "yugd@live.cn",
        license = "GPL",
        keywords = ("Guodong", "high-throughput"),
        platforms = "Independant",
        url = "",
      )

