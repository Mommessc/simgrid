/* Copyright (c) 2016-2020. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SIMGRID_PLUGINS_TEMPERATURE_H_
#define SIMGRID_PLUGINS_TEMPERATURE_H_

#include <xbt/base.h>
#include <simgrid/forward.h>

SG_BEGIN_DECL

XBT_PUBLIC void sg_host_temperature_plugin_init();
XBT_PUBLIC double sg_host_get_aggregated_consumed_energy(const_sg_host_t host);
XBT_PUBLIC double sg_host_get_host_temperature(const_sg_host_t host);
XBT_PUBLIC double sg_host_get_air_temperature(const_sg_host_t host);
XBT_PUBLIC double sg_host_get_outside_temperature(const_sg_host_t host);
XBT_PUBLIC void sg_host_set_outside_temperature(const_sg_host_t host, double temperature);
XBT_PUBLIC void sg_host_update_temperatures(const_sg_host_t host);
XBT_PUBLIC void sg_host_init_master_slaves(const_sg_host_t host);



#define MSG_host_temperature_plugin_init() sg_host_temperature_plugin_init()
#define MSG_host_get_aggregated_consumed_energy(host) sg_host_get_aggregated_consumed_energy(host)
#define MSG_host_get_host_temperature(host) sg_host_get_host_temperature(host)
#define MSG_host_get_air_temperature(host) sg_host_get_air_temperature(host)
#define MSG_host_get_outside_temperature(host) sg_host_get_outside_temperature(host)
#define MSG_host_set_outside_temperature(host, temperature) sg_host_set_outside_temperature(host, temperature)
#define MSG_host_update_temperatures(host) sg_host_update_temperatures(host)
#define MSG_host_init_master_slaves(host) sg_host_init_master_slaves(host)

SG_END_DECL

#endif
