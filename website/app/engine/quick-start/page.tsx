"use client";

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { motion } from "framer-motion";

export default function QuickStartPage() {
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
            Quick Start Guide
          </span>
        </h1>
        <p className="text-xl text-muted-foreground">
          Get started with the fuel blend optimization engine in minutes.
        </p>
      </motion.div>

      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.2 }}
        className="space-y-8"
      >
        <motion.section
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 0.4 }}
        >
          <h2 className="mb-4 text-3xl font-bold">Installation</h2>
          <Card>
            <CardContent className="pt-6">
              <pre className="overflow-x-auto rounded-lg bg-muted p-4 text-sm">
                <code>{`# Install dependencies
pip install -r requirements.txt`}</code>
              </pre>
            </CardContent>
          </Card>
        </motion.section>

        <motion.section
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 0.6 }}
        >
          <h2 className="mb-4 text-3xl font-bold">Basic Usage</h2>
          <Card>
            <CardHeader>
              <CardTitle>Command Line Interface</CardTitle>
            </CardHeader>
            <CardContent className="space-y-4">
              <div>
                <p className="mb-2 font-semibold">Gasoline Optimization</p>
                <pre className="overflow-x-auto rounded-lg bg-muted p-4 text-sm">
                  <code>{`python engine.py \\
  --components data/components.csv \\
  --fuel-type gasoline \\
  --K 5 \\
  --top 15 \\
  --seed 42`}</code>
                </pre>
              </div>
              <div>
                <p className="mb-2 font-semibold">Diesel Optimization</p>
                <pre className="overflow-x-auto rounded-lg bg-muted p-4 text-sm">
                  <code>{`python engine.py \\
  --components data/components.csv \\
  --fuel-type diesel \\
  --K 4 \\
  --top 15 \\
  --seed 123`}</code>
                </pre>
              </div>
              <div>
                <p className="mb-2 font-semibold">With Novel Components</p>
                <pre className="overflow-x-auto rounded-lg bg-muted p-4 text-sm">
                  <code>{`python engine.py \\
  --components data/components.csv \\
  --novel data/components_novel.csv \\
  --allow-novel 1 \\
  --max-add-novel 0.1 \\
  --fuel-type gasoline \\
  --K 5 \\
  --top 10`}</code>
                </pre>
              </div>
            </CardContent>
          </Card>
        </motion.section>

        <motion.section
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 0.8 }}
        >
          <h2 className="mb-4 text-3xl font-bold">Key Parameters</h2>
          <div className="grid gap-4 md:grid-cols-2">
            <Card>
              <CardHeader>
                <CardTitle>--components</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  Path to component CSV file (required). Must include columns: id, name, density_g_ml,
                  LHV_MJ_kg, RON, MON, cetane, PMI, TSI, cost_eur_L, max_vol_frac, is_novel.
                </p>
              </CardContent>
            </Card>
            <Card>
              <CardHeader>
                <CardTitle>--fuel-type</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  Fuel type: gasoline, diesel, or jet. Determines which specifications are enforced.
                </p>
              </CardContent>
            </Card>
            <Card>
              <CardHeader>
                <CardTitle>--K</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  Maximum number of components in the blend (default: 5).
                </p>
              </CardContent>
            </Card>
            <Card>
              <CardHeader>
                <CardTitle>--top</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  Number of top results to return (default: 15).
                </p>
              </CardContent>
            </Card>
            <Card>
              <CardHeader>
                <CardTitle>--seed</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  Random seed for reproducibility (default: 42).
                </p>
              </CardContent>
            </Card>
            <Card>
              <CardHeader>
                <CardTitle>--allow-novel</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  Allow novel components (0 or 1). Requires --novel file to be specified.
                </p>
              </CardContent>
            </Card>
          </div>
        </motion.section>

        <motion.section
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 1.0 }}
        >
          <h2 className="mb-4 text-3xl font-bold">Output</h2>
          <Card>
            <CardContent className="pt-6">
              <p className="mb-4 text-muted-foreground">
                The engine generates output files in <code className="rounded bg-muted px-1">runs/&lt;timestamp&gt;_&lt;seed&gt;/</code>:
              </p>
              <ul className="space-y-2 text-muted-foreground">
                <li className="flex items-start gap-2">
                  <span className="mt-1 text-primary">•</span>
                  <span><strong className="text-foreground">Top-N.csv:</strong> Optimized blends with properties, component fractions, lab recipes, and cost analysis</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="mt-1 text-primary">•</span>
                  <span><strong className="text-foreground">RUN.json:</strong> Complete run metadata including parameters and configuration</span>
                </li>
              </ul>
            </CardContent>
          </Card>
        </motion.section>
      </motion.div>
    </div>
  );
}

