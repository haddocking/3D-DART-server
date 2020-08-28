FROM ubuntu:18.04

MAINTAINER BonvinLab <software.csb@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive 

RUN apt-get update && apt-get install -y apache2 python2.7 python2.7-numpy python-pip vim gcc-multilib
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

ENV APACHE_RUN_USER  www-data
ENV APACHE_RUN_GROUP www-data
ENV APACHE_LOG_DIR   /var/log/apache2
ENV APACHE_PID_FILE  /var/run/apache2/apache2.pid
ENV APACHE_RUN_DIR   /var/run/apache2
ENV APACHE_LOCK_DIR  /var/lock/apache2
ENV APACHE_LOG_DIR   /var/log/apache2

RUN mkdir -p $APACHE_RUN_DIR
RUN mkdir -p $APACHE_LOCK_DIR
RUN mkdir -p $APACHE_LOG_DIR

COPY apache2 /etc/apache2
RUN ln -s /etc/apache2/mods-available/cgi.load /etc/apache2/mods-enabled/cgi.load

RUN mkdir -p /var/www/html/3DDART
COPY --chown=www-data:www-data 3d-dart /var/www/html/3DDART/
RUN mkdir -p /var/www/html/3DDART/error && chown -R www-data:www-data /var/www/html/3DDART/error
RUN ln -s /var/www/html/3DDART/server/server-tmp/software /var/www/html/3DDART/software
RUN cp /var/www/html/3DDART/index.html /var/www/html/index.html

ENV X3DNA=/var/www/html/3DDART/server/server-tmp/software/X3DNA-linux
ENV PATH=$PATH:/var/www/html/3DDART/server/server-tmp/software/X3DNA-linux/bin:/var/www/html/3DDART/server

EXPOSE 80

CMD ["/usr/sbin/apache2", "-D", "FOREGROUND"]