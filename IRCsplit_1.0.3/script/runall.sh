#!/bin/bash
for inf in *.gjf
do
g09 < ${inf} |tee ${inf//gjf/out}
done
