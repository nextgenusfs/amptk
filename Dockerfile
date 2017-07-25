FROM amptk

WORKDIR /amptk

COPY usearch9.2.64_i86linux32 /amptk/usearch9

RUN chmod +x /amptk/usearch9

WORKDIR /work
