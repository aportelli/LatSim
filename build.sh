#!/bin/bash

set -e

PREFIX=`cat Makefile | grep '^prefix =' | awk '{print $3}'`
case $1 in
    '')
        echo '-- building...'
        make uninstall 1>/dev/null
        make -j8
        echo '-- installing...'
        make install 1>/dev/null;;
    'clean')
        echo '-- cleaning...'
        make -j8 clean;;
    *)
        echo 'error: unknown action' 1>&2
        exit 1;;
esac
