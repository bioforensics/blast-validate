#! /usr/bin/perl

use strict;

while ( <> )
{
    s/^@/>/;
    print;
    print scalar <>;
    <>;
    <>;
}
