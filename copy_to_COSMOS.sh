#! bash

for exe in bam-quant bam-coverage bw-compare; do
    cp "target/x86_64-unknown-linux-musl/release/$exe" \
       "/home/med-sal/sens05_home/bin/"
done
