#!/bin/sh
/usr/bin/x11vnc -xkb -noxrecord -noxfixes -noxdamage -display :0 -repeat -nopw -wait 5 -permitfiletransfer -tightfilexfer -auth /var/run/slim.auth
#/usr/bin/supervisorctl restart slim


