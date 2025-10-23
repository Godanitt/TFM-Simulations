#!/bin/bash
set -e

# ======================================
# CONFIGURACIÓN GENERAL
# ======================================
INSTALL_DIR=$HOME
NPROC=$(nproc)

echo "====> Instalando dependencias básicas del sistema..."
sudo apt-get update
sudo apt-get install -y \
  build-essential gfortran git wget curl cmake cmake-curses-gui \
  libssl-dev libx11-dev libxpm-dev libxft-dev libxext-dev \
  libxmu-dev libxi-dev libxrender-dev libxrandr-dev libxt-dev \
  libglu1-mesa-dev freeglut3-dev mesa-common-dev libglew-dev \
  libfftw3-dev libcfitsio-dev libpng-dev libjpeg-dev \
  qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools \
  libexpat1-dev libxml2-dev \
  python3 python3-pip \
  unzip

echo "====> Todas las dependencias del sistema instaladas."

# ======================================
# INSTALAR CMAKE (si se desea más nuevo)
# ======================================
sudo apt-get install -y cmake

# ======================================
# INSTALAR ROOT
# ======================================
cd $INSTALL_DIR
echo "====> Descargando ROOT..."
git clone --branch latest-stable https://github.com/root-project/root.git
mkdir -p root-build root-install
cd root-build

echo "====> Configurando ROOT..."
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/root-install \
      -DCMAKE_CXX_STANDARD=17 \
      -Dgnuinstall=ON \
      ../root

echo "====> Compilando ROOT..."
make -j$NPROC
make install
echo "====> ROOT instalado en $INSTALL_DIR/root-install"

# ======================================
# INSTALAR GEANT4
# ======================================
cd $INSTALL_DIR
echo "====> Descargando Geant4..."
wget https://gitlab.cern.ch/geant4/geant4/-/archive/v11.3.2/geant4-v11.3.2.tar.gz
tar -xzf geant4-v11.3.2.tar.gz
mkdir -p geant4-build geant4-install
cd geant4-build

echo "====> Configurando Geant4..."
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/geant4-install \
      -DGEANT4_INSTALL_DATA=ON \
      -DGEANT4_USE_GDML=ON \
      -DGEANT4_USE_OPENGL_X11=ON \
      -DGEANT4_USE_QT=ON \
      ../geant4-v11.3.2

echo "====> Compilando Geant4..."
make -j$NPROC
make install
echo "====> Geant4 instalado en $INSTALL_DIR/geant4-install"

# ======================================
# INSTALAR GARFIELD++
# ======================================
cd $INSTALL_DIR
echo "====> Descargando Garfield++..."
git clone https://gitlab.cern.ch/garfield/garfieldpp.git garfield
mkdir -p garfield-build garfield-install
cd garfield-build

echo "====> Configurando Garfield++..."
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/garfield-install \
      -DGeant4_DIR=$INSTALL_DIR/geant4-install/lib/cmake/Geant4 \
      -DROOT_DIR=$INSTALL_DIR/root-install/cmake \
      ../garfield

echo "====> Compilando Garfield++..."
make -j$NPROC
make install
echo "====> Garfield++ instalado en $INSTALL_DIR/garfield-install"

# ======================================
# CONFIGURAR VARIABLES DE ENTORNO
# ======================================
echo "====> Añadiendo variables de entorno al ~/.bashrc ..."
{
  echo ""
  echo "# ==== Entorno ROOT, Geant4, Garfield++ ===="
  echo "source $INSTALL_DIR/root-install/bin/thisroot.sh"
  echo "source $INSTALL_DIR/geant4-install/bin/geant4.sh"
  echo "export LD_LIBRARY_PATH=$INSTALL_DIR/garfield-install/lib:\$LD_LIBRARY_PATH"
  echo "export PATH=$INSTALL_DIR/garfield-install/bin:\$PATH"
} >> ~/.bashrc

echo "====> Instalación completada con éxito."
echo "Reinicia la terminal o ejecuta: source ~/.bashrc"
