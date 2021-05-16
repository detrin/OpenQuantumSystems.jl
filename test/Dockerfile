FROM python:3.9

ENV PYTHONUNBUFFERED=1

RUN python -m ensurepip
RUN python -m pip install --no-cache --upgrade pip setuptools

RUN python -m pip install quantarhei

WORKDIR /root/app

CMD ["python", "hamiltonian.py"]