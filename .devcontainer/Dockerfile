FROM python:3.9

# Set environment variables to non-interactive (this prevents some prompts)
ENV DEBIAN_FRONTEND=non-interactive
ENV PYTHONUNBUFFERED=true

RUN apt-get update && \
    apt-get install -y build-essential libffi-dev && \
	apt-get remove && \
	apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip
RUN pip install poetry

RUN poetry config virtualenvs.in-project true