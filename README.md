# Extrapolation

Extrapolating the final results between signal lifetimes. 

### Clone repository

```bash
git clone https://github.com/apmorris/Extrapolation.git
```

### Set up (lxplus)

```bash
setupATLAS
lsetup git
lsetup "root 6.10.04-x86_64-slc6-gcc62-opt"
# Make sure local copy is up to date
git pull origin master
```

### Run limit finder

```bash
cd FindLimit/
make clean
make 
./FindLimit -A <nA> -B <nB> -C <nC> -D <nD> -w <sA> -x <sB> -y <sC> -z <sD> (-a)
```
_NB:_ systematic errors are currently hardcoded into LimitCommonCode/run_ABCD.cxx

### Run extrapolation

```bash
cd ExtrapLimitFinder/
make clean
make
./ExtrapLimitFinder
```

### Commit any changes

```bash
git status
git add <files>
git commit -m "Descriptive commit message"
git push origin master
```

Code inherited from Gordon Watts, modified to work on lxplus.