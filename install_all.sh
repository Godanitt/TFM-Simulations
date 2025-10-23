#!/bin/bash
set -e
NPROC=$(nproc)
INSTALL_DIR=$HOME

echo "=============================="
echo " Instalador ROOT + Geant4 + Garfield++"
echo " Ubuntu 20.04 - configurado automáticamente"
echo "=============================="

# --- Dependencias del sistema ---
echo "[1/8] Instalando dependencias del sistema..."
sudo apt-get update
sudo apt-get install -y \
  build-essential gfortran git wget curl cmake cmake-curses-gui \
  libssl-dev libx11-dev libxpm-dev libxft-dev libxext-dev \
  libxmu-dev libxi-dev libxrender-dev libxrandr-dev libxt-dev \
  libglu1-mesa-dev freeglut3-dev mesa-common-dev libglew-dev \
  libfftw3-dev libcfitsio-dev libpng-dev libjpeg-dev \
  qt5-default qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools \
  libexpat1-dev libxml2-dev python3 python3-pip unzip

# --- Instalar ROOT ---
echo "[2/8] Descargando ROOT..."
cd $INSTALL_DIR
git clone --branch v6-30-02 https://github.com/root-project/root.git
mkdir -p root-build root-install
cd root-build

echo "[3/8] Configurando ROOT..."
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/root-install \
      -DCMAKE_CXX_STANDARD=17 \
      -Dgnuinstall=ON \
      -Dx11=ON \
      ../root

echo "[4/8] Compilando ROOT (esto tardará varios minutos)..."
make -j$NPROC
make install

# --- Instalar Geant4 ---
echo "[5/8] Descargando Geant4..."
cd $INSTALL_DIR
wget https://gitlab.cern.ch/geant4/geant4/-/archive/v11.2.2/geant4-v11.2.2.tar.gz
tar -xzf geant4-v11.2.2.tar.gz
mkdir -p geant4-build geant4-install
cd geant4-build

echo "[6/8] Configurando Geant4..."
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/geant4-install \
      -DGEANT4_INSTALL_DATA=ON \
      -DGEANT4_USE_GDML=ON \
      -DGEANT4_USE_OPENGL_X11=ON \
      -DGEANT4_USE_QT=ON \
      ../geant4-v11.2.2

echo "[7/8] Compilando e instalando Geant4..."
make -j$NPROC
make install

# --- Instalar Garfield++ ---
echo "[8/8] Descargando y compilando Garfield++..."
cd $INSTALL_DIR
git clone https://gitlab.cern.ch/garfield/garfieldpp.git garfield
mkdir -p garfield-build garfield-install
cd garfield-build
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/garfield-install \
      -DGeant4_DIR=$INSTALL_DIR/geant4-install/lib/cmake/Geant4 \
      -DROOT_DIR=$INSTALL_DIR/root-install/cmake \
      ../garfield
make -j$NPROC
make install

# --- Configurar entorno ---
echo "Configurando entorno permanente en ~/.bashrc ..."
{
  echo ""
  echo "# ==== Entorno ROOT, Geant4, Garfield++ ===="
  echo "source $INSTALL_DIR/root-install/bin/thisroot.sh"
  echo "source $INSTALL_DIR/geant4-install/bin/geant4.sh"
  echo "export LD_LIBRARY_PATH=$INSTALL_DIR/garfield-install/lib:\$LD_LIBRARY_PATH"
  echo "export PATH=$INSTALL_DIR/garfield-install/bin:\$PATH"
} >> ~/.bashrc

echo "=============================="
echo " Instalación completada 🎉"
echo " Ejecuta ahora: source ~/.bashrc"
echo "=============================="
