"use client";

import { PieChart, Pie, Cell, ResponsiveContainer, Legend, Tooltip } from "recharts";

const COLORS = [
  "oklch(0.75 0.22 50)", // Yellow
  "oklch(0.70 0.20 40)", // Amber
  "oklch(0.65 0.20 35)", // Orange
  "oklch(0.60 0.18 25)", // Red-orange
  "oklch(0.35 0.10 15)", // Deep red
];

interface BlendCompositionProps {
  components: Array<{ id: string; name: string; vol_frac: number }>;
}

export function BlendComposition({ components }: BlendCompositionProps) {
  const data = components
    .filter((c) => c.vol_frac > 0)
    .map((c) => ({
      name: c.name || c.id,
      value: c.vol_frac * 100,
    }));

  return (
    <ResponsiveContainer width="100%" height={300}>
      <PieChart>
        <Pie
          data={data}
          cx="50%"
          cy="50%"
          labelLine={false}
          label={({ name, percent }) => `${name}: ${((percent ?? 0) * 100).toFixed(0)}%`}
          outerRadius={80}
          fill="#8884d8"
          dataKey="value"
        >
          {data.map((entry, index) => (
            <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
          ))}
        </Pie>
        <Tooltip
          contentStyle={{
            backgroundColor: "oklch(0.10 0 0)",
            border: "1px solid oklch(0.65 0.20 35)",
            borderRadius: "0.5rem",
          }}
        />
        <Legend />
      </PieChart>
    </ResponsiveContainer>
  );
}

