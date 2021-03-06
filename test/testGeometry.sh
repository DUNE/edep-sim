#!/bin/bash

#!/bin/sh

OUTPUT=testGeometry.root

if [ -f ${OUTPUT} ]; then
    rm ${OUTPUT}
fi

MACRO=$(mktemp -t detsim.XXXXXXXXXX) || exit 1
cat > ${MACRO} <<EOF
/gps/particle mu+
/gps/energy 100 MeV
/gps/position 0.0 0.0 -15.0 cm
/gps/pos/type Volume
/gps/ang/type iso
EOF

edep-sim -o ${OUTPUT} -C -e 1 -u ${MACRO}

