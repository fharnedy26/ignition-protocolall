"use client";

import Link from "next/link";
import { ArrowRight, BookOpen, Zap, Target } from "lucide-react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { motion } from "framer-motion";

export default function ProjectPage() {
  return (
    <div className="container mx-auto px-4 py-12">
      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6 }}
        className="mb-12"
      >
        <h1 className="mb-4 text-5xl font-bold">
          <span className="bg-gradient-to-r from-primary to-accent bg-clip-text text-transparent">
            Ignition Protocol
          </span>
        </h1>
        <p className="text-xl text-muted-foreground">
          An Evolutionary Engine for Molecular Combustion
        </p>
        <p className="mt-2 text-sm text-muted-foreground">
          Stripe Young Scientist & Technology Exhibition 2026 • by Finn Harnedy
        </p>
      </motion.div>

      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.2 }}
        className="mb-12"
      >
        <div className="mb-6">
          <h2 className="mb-4 text-3xl font-bold">Project Overview</h2>
          <p className="text-muted-foreground">
            Ignition Protocol is a computational framework for discovering, filtering, and operationalising 
            synthetic fuel candidates. The system combines molecular generation with blend optimization to 
            create an end-to-end fuel design pipeline.
          </p>
        </div>
        <motion.div
          whileHover={{ scale: 1.02, transition: { duration: 0.2 } }}
        >
          <Card className="border-primary/50 bg-primary/5">
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <BookOpen className="h-5 w-5 text-primary" />
                Full Project Book Available
              </CardTitle>
            </CardHeader>
            <CardContent>
              <p className="text-muted-foreground">
                The complete Ignition Protocol project book with detailed methodology, results, and analysis 
                is available at the exhibition stand. Visit us to learn more!
              </p>
            </CardContent>
          </Card>
        </motion.div>
      </motion.div>

      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.4 }}
        className="mb-12"
      >
        <h2 className="mb-6 text-3xl font-bold">Key Highlights</h2>
        <div className="grid gap-6 md:grid-cols-2">
          <Link href="/project/highlights">
            <motion.div
              whileHover={{ scale: 1.03, y: -5 }}
              transition={{ duration: 0.2 }}
            >
              <Card className="group h-full transition-all hover:border-primary hover:shadow-[0_0_20px_oklch(0.65_0.20_35/0.3)]">
                <CardHeader>
                  <div className="mb-2 flex h-12 w-12 items-center justify-center rounded-lg bg-primary/10 text-primary">
                    <Zap className="h-6 w-6" />
                  </div>
                  <CardTitle>Key Highlights</CardTitle>
                </CardHeader>
                <CardContent>
                  <p className="text-muted-foreground">
                    Explore the abstract, key methodology points, notable results, and future directions.
                  </p>
                  <div className="mt-4 flex items-center gap-2 text-primary group-hover:gap-4 transition-all">
                    <span className="text-sm font-medium">View highlights</span>
                    <ArrowRight className="h-4 w-4" />
                  </div>
                </CardContent>
              </Card>
            </motion.div>
          </Link>
          <Link href="/project/methodology">
            <motion.div
              whileHover={{ scale: 1.03, y: -5 }}
              transition={{ duration: 0.2 }}
            >
              <Card className="group h-full transition-all hover:border-primary hover:shadow-[0_0_20px_oklch(0.65_0.20_35/0.3)]">
                <CardHeader>
                  <div className="mb-2 flex h-12 w-12 items-center justify-center rounded-lg bg-primary/10 text-primary">
                    <Target className="h-6 w-6" />
                  </div>
                  <CardTitle>Methodology</CardTitle>
                </CardHeader>
                <CardContent>
                  <p className="text-muted-foreground">
                    Learn about the evolutionary algorithms, fitness functions, and optimization approach.
                  </p>
                  <div className="mt-4 flex items-center gap-2 text-primary group-hover:gap-4 transition-all">
                    <span className="text-sm font-medium">View methodology</span>
                    <ArrowRight className="h-4 w-4" />
                  </div>
                </CardContent>
              </Card>
            </motion.div>
          </Link>
        </div>
      </motion.div>

      <motion.section
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.6 }}
      >
        <h2 className="mb-4 text-3xl font-bold">System Architecture</h2>
        <Card>
          <CardContent className="pt-6">
            <p className="mb-4 text-muted-foreground">
              Ignition Protocol is a unified pipeline with two coupled engines:
            </p>
            <div className="space-y-4">
              <motion.div
                initial={{ opacity: 0, x: -20 }}
                animate={{ opacity: 1, x: 0 }}
                transition={{ duration: 0.5, delay: 0.8 }}
                className="rounded-lg border border-border bg-card p-4"
              >
                <h3 className="mb-2 font-semibold text-primary">1. Molecular Generation Module</h3>
                <p className="text-muted-foreground">
                  Evolves candidate molecules from fragment libraries using RDKit-based validity checks 
                  and a multi-objective fitness score. Produces candidate SMILES strings with practical 
                  metadata including estimated fuel properties, synthesis feasibility, and storage stability proxies.
                </p>
              </motion.div>
              <motion.div
                initial={{ opacity: 0, x: -20 }}
                animate={{ opacity: 1, x: 0 }}
                transition={{ duration: 0.5, delay: 1.0 }}
                className="rounded-lg border border-border bg-card p-4"
              >
                <h3 className="mb-2 font-semibold text-primary">2. Blend Optimization Engine</h3>
                <p className="text-muted-foreground">
                  Deterministic optimization that selects and proportions real components (including optional 
                  novel candidates) into lab-ready fuel blends. Enforces constraints such as volatility limits 
                  (T10/T90, RVP), minimum octane/cetane requirements, density bounds, and component caps.
                </p>
              </motion.div>
            </div>
          </CardContent>
        </Card>
      </motion.section>
    </div>
  );
}

