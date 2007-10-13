
#!/bin/bash
# vim: set ft=sh sw=4 sts=4 et :

echo ">>> test ${2:-${1}}"
${@} quick mc
