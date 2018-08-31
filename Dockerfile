FROM gradle:jdk8 as builder
COPY --chown=gradle:gradle ./src /home/gradle/src
COPY --chown=gradle:gradle ./lib /home/gradle/lib
COPY --chown=gradle:gradle build.gradle settings.gradle /home/gradle/
WORKDIR /home/gradle/
RUN gradle build
RUN ls build/distributions

FROM openjdk:8-jre-slim
EXPOSE 9090
COPY --from=builder /home/gradle/build/distributions/vassar_server-1.0.tar /app/
COPY ./problems /app/problems
COPY ./resources  /app/resources
COPY ./templates /app/templates
WORKDIR /app/
RUN tar -xvf vassar_server-1.0.tar
RUN ls vassar_server-1.0/bin
CMD vassar_server-1.0/bin/vassar_server