/*
zon wil print timestamp for sun rise and set
Copyright (C) 2022 Copyright Michael Welter

This file is part "zon"

zon is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

zon is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with zon.  If not, see <http://www.gnu.org/licenses/>.

*/

/* The options we understand. */
static struct argp_option options[] =
{
    {0,0,0,0, "Options to select what to display" },
    {"rise",     'r', 0,      0,  "Produce the next start/rise time relative to current or provided time" },
    {"set",      's', 0,      0,  "Produce the next end/set time relative to current or provided time" },
    {0,0,0,0, "" },
    {"mid",      'm', 0,      0,  "Produce time exactly between the next rise and set times, i.e. deep midnight or high noon" },
    {"current",  'c', 0,      0,  "Whether sun is up \"+\" or down \"-\". Default if no display type selected." },
    {0,0,0,0, "Options to specify when and where on earth" },
    {"location", 'l', "+DDMM+DDDMM", 0,  "or +DDMMSS+DDDMMSS or degrees,degrees (with N,S,+ or - signs) Calculate for latitude (+N/-S) and longitude (+E-W) in Degrees, Minutes and Seconds. Overrides configuration files /etc/zon.conf and ~/.config/zon.conf" },
    {"date",     'd', "iso-time", 0,  "YYYY-MM-DDTHH:MM+ZZ[:]zz Calculate for specified iso-formatted time. Defaults to current system time. Specify \"date <date syntax>\" to parse by invoking the date command" },
    {"file",  'f', "file",      0,  "like --date; once for each line of FILE. \"-\" for stdin" },
    {0,0,0,0, "Options to format the output. Defaults to iso-format" },
    {"at",       '@', 0,      0,  "Format output as date usable by the at command (in UTC), HH:MM YYYY-MM-DD" },
    {"systemd",  'y', 0,      0,  "Format output as required for systemd-run,  YYYY-MM-DD HH:MM UTC" },
    {"systemdfile",  'Y', 0,      0,  "Format output as required for systemd.timer, [Timer]\\nOnCalendar=YYYY-MM-DD HH:MM UTC" },
    {0,0,0,0, "" },
    {"verbose",  'v', 0,      0,  "Produce a label or if repeated give all base and calculated data, including date and location" },
    {"output",  'o', "filename",      0,  "Output to file, replacing any previous content" },
    {0,0,0,0, "Options to select the kind of twighlight" },
    {"sun",       0,  0,      0,  "Default: Produce start, ending and duration of visibility of top of sun above horizon, i.e. sunrise and sunset. Both atmospheric refraction (-35/60 degree) and rim of the apparent size of the solar disk are accounted for." },
    {0,0,0,0, "" },
    {"civil",     1,  0,      0,  "Produce data about civil twighlight, starting when centre of sun is 6 degrees below horizon" },
    {0,0,0,0, "" },
    {"nautical",  2,  0,      0,  "Produce data about nautical twighlight, starting when centre of sun is 12 degrees below horizon" },
    {0,0,0,0, "" },
    {"astronomical",3,0,      0,  "Produce data about astronomical twighlight, starting when centre of sun is 18 degrees below horizon" },
    {0,0,0,0, "" },
    {"angle"     ,5,  "degrees",      0,  "Specify your own rise/set angle of the centre of the sun to the horizon" },
    {"rim"       ,4,  0,      0,  "Specify to compensate angle for the upper rim of the sun (i.e. the radius of the apparent solar disk). Like at sun rise/set. Use after --angle" },
    {"phase"     ,6,  0,      0,  "Produce data about phase of the moon" },
    { 0 }
};
