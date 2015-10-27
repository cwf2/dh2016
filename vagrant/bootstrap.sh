# Install prerequisites on Vagrant VM
# using ubuntu/trusty64

sudo -u vagrant cp /vagrant/vagrant/screenrc /home/vagrant/.screenrc
sudo -u vagrant cp /vagrant/vagrant/vimrc /home/vagrant/.vimrc

apt-get install -y \
  git \
  htop \
  libxml2-dev

# set up R

sudo -u vagrant gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
sudo -u vagrant bash -c "gpg -a --export E084DAB9 | sudo apt-key add -"

echo http://cran.rstudio.com/bin/linux/ubuntu trusty/ >> /etc/apt/source.list
  apt-get update
  apt-get install -y r-base-dev

sudo -u vagrant bash -c "R --vanilla < /vagrant/vagrant/install.R"

