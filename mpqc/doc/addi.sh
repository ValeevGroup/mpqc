#!/bin/sh

echo -n ' -I '`echo $* | sed "s/ / -I /"`
