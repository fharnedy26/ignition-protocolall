# Debian 12 Installation Guide for MacBook Pro 2010 17-inch (i7)

## Your MacBook Pro 2010 17-inch Specifications

**Hardware:**
- **CPU**: Intel Core i7 (64-bit, use amd64 ISO)
- **WiFi**: Broadcom BCM4322 (needs special drivers)
- **Graphics**: NVIDIA GeForce 330M + Intel HD Graphics (dual GPU)
- **Audio**: NVIDIA MCP89 (works with snd-hda-intel)
- **Bluetooth**: Broadcom BCM2046 (works with btusb)
- **Boot**: EFI (64-bit)

**Known Issues:**
- WiFi may need proprietary driver (broadcom-sta-dkms) for best results
- WiFi interface may be "dormant" after boot (sleep/wake workaround)
- Graphics switching may need configuration

## Pre-Installation Checklist

### 1. Download Debian 12 ISO

**For MacBook Pro 2010 17-inch:**

**CRITICAL: Use the firmware ISO with non-free drivers!**

- **Download**: Debian 12 netinst ISO **with firmware** (includes WiFi drivers)
  - URL: https://cdimage.debian.org/cdimage/unofficial/non-free/cd-including-firmware/
  - Choose: `firmware-12.x.x-amd64-netinst.iso`
  - **This is essential** - your BCM4322 WiFi card needs these drivers

**Alternative (if firmware ISO doesn't work):**
- Standard netinst ISO: https://www.debian.org/CD/netinst/
- Choose: `debian-12.x.x-amd64-netinst.iso`
- You'll need to install WiFi drivers manually after installation

**Note**: The "netinst" (network install) version is smaller (~500MB) and downloads packages during installation. Perfect for MacBooks.

### 2. Create Bootable USB Stick

**On macOS (before installing Debian):**

#### Method 1: Using `dd` command (Terminal)

```bash
# 1. Find your USB stick device name
diskutil list
# Look for your USB stick (usually /dev/disk2 or /dev/disk3)
# Note: Use /dev/rdiskX (with 'r') for faster writes

# 2. Unmount the USB stick (replace /dev/disk2 with your device)
diskutil unmountDisk /dev/disk2

# 3. Write the ISO to USB (replace paths as needed)
sudo dd if=/path/to/firmware-12.x.x-amd64-netinst.iso of=/dev/rdisk2 bs=1m
# Wait for it to complete (can take 5-10 minutes)
# You'll see "bytes transferred" when done

# 4. Eject the USB stick
diskutil eject /dev/disk2
```

#### Method 2: Use Etcher (easier, GUI tool - RECOMMENDED)
- Download: https://etcher.balena.io/
- Select ISO file
- Select USB stick
- Flash
- Verify (Etcher does this automatically)

#### Method 3: Using Disk Utility (macOS built-in)
1. Open Disk Utility
2. Select your USB stick
3. Click "Erase" → Format as "MS-DOS (FAT32)" or "ExFAT"
4. **Then use `dd` command above** (Disk Utility alone won't create bootable USB from ISO)

### 3. Verify USB Drive After Flashing

**What should be on the USB drive:**

After flashing, the USB drive should contain:
- **EFI folder** (for EFI boot)
- **boot folder** (bootloader files)
- **dists folder** (package lists)
- **pool folder** (packages)
- **debian folder** (installer files)
- Various configuration files (.disk, README, etc.)

**To verify on macOS:**

```bash
# 1. Check if USB is mounted
diskutil list

# 2. Check the partition (USB should show as ISO 9660 or similar)
diskutil info /dev/disk2s1  # Replace with your USB device

# 3. Try to mount and view contents
# Note: macOS may not mount it properly - this is normal
# The important thing is it boots on the MacBook
```

**Important Notes:**
- macOS may show the USB as "unreadable" or "not mounted" - **this is normal**
- The ISO creates a hybrid filesystem that macOS doesn't fully understand
- **What matters is if it boots on the MacBook**, not if macOS can read it
- If the MacBook doesn't see it, the flash may have failed

### 3. MacBook-Specific Preparations

**Before installation:**

1. **Backup your data** (if dual-booting or replacing macOS)
   - Your 2010 MacBook Pro doesn't have FileVault, but still backup important data

2. **MacBook Pro 2010 specific notes**:
   - WiFi: Broadcom BCM4322 - will need proprietary driver after installation
   - Graphics: NVIDIA GeForce 330M - may need NVIDIA drivers (optional, Intel works fine)
   - Audio: Should work out of the box
   - Bluetooth: Should work out of the box
   - Suspend: Works out of the box

3. **Use firmware ISO** (recommended):
   - Download the firmware ISO from the link above
   - This includes WiFi drivers that will help during installation

4. **Boot considerations**:
   - Your 2010 MacBook uses EFI boot (64-bit)
   - Hold Option key during boot to select USB
   - No Secure Boot (not available on 2010 models)

## Installation Steps

### Step 1: Boot from USB

**If USB drive is not detected:**

1. **Verify the flash was successful:**
   ```bash
   # On macOS, check the USB was written to
   diskutil list
   # The USB should show a partition
   ```

2. **Try different USB ports:**
   - Use a USB 2.0 port if available (not USB 3.0)
   - Try the USB port on the opposite side
   - Some MacBooks are picky about which port works

3. **Boot process:**
   - **Shut down your MacBook completely**
   - **Insert the USB stick**
   - **Power on and immediately hold Option (Alt) key**
   - **Keep holding until you see boot menu** (may take 10-15 seconds)
   - **Look for**: "EFI Boot", "Windows", or the USB stick name
   - **If nothing appears**, try:
     - Hold Option longer (up to 30 seconds)
     - Try a different USB stick
     - Try re-flashing the ISO

4. **If USB still doesn't appear:**
   - The flash may have failed - try re-flashing
   - Use a different USB stick (some don't work well)
   - Try Etcher instead of `dd` command
   - Make sure you used `/dev/rdiskX` not `/dev/diskX` (the 'r' is important)

5. **Once USB appears in boot menu:**
   - **Select the USB stick** (usually labeled "EFI Boot" or "Windows")
   - **Select "Install" or "Graphical Install"**

### Step 2: Installation Process

**Key choices during installation:**

1. **Language**: Choose your language
2. **Location**: Select your country/region
3. **Keyboard**: 
   - For MacBook: Select "Macintosh" or "Apple" keyboard layout
   - Or: "English (US)" if you prefer standard layout
4. **Network Configuration** (IMPORTANT FOR MACBOOK PRO 2010):
   - **WiFi will NOT work during installation** (drivers not installed yet)
   - **You'll see "Choose network" with no options** - this is normal!
   - **Select "Do not configure the network at this time"** or similar option
   - **OR** if you have Ethernet cable, connect it and use that
   - You'll configure WiFi AFTER installation is complete
   - Installation can continue without network (packages will be downloaded later)
5. **Partitioning**:
   - **Option A (Recommended)**: "Guided - use entire disk"
     - This will erase everything on the disk
   - **Option B**: "Guided - use entire disk and set up LVM"
     - More flexible for future changes
   - **Option C**: Manual (advanced users only)
6. **Software Selection**:
   - ✅ **Check "Debian desktop environment"** (this is the windowing system; needed if you want to use a graphical web browser to search the website)
   - ✅ **Check "SSH server"** (for remote access)
   - ✅ **Check "standard system utilities"**
   - You can install more later
7. **Boot Loader**: Install GRUB to the master boot record (default is fine)

### Step 3: Post-Installation (First Boot)

**After installation completes and system reboots:**

1. **Log in** with the root password you set
2. **Create a user account** (if you didn't during installation):
   ```bash
   adduser yourusername
   usermod -aG sudo yourusername
   ```

3. **Update the system**:
   ```bash
   apt update
   apt upgrade -y
   ```

## MacBook-Specific Setup

### WiFi Configuration (Broadcom BCM4322 - MacBook Pro 2010)

**Your MacBook Pro 2010 has a Broadcom BCM4322 WiFi card. You have two driver options:**

#### Option 1: Proprietary `wl` Driver (RECOMMENDED - Better Performance)

```bash
# 1. Add non-free repositories
sudo nano /etc/apt/sources.list
# Add "non-free" and "contrib" to each line, e.g.:
# deb http://deb.debian.org/debian bookworm main non-free contrib

# 2. Update package list
sudo apt update

# 3. Install the proprietary Broadcom driver
sudo apt install broadcom-sta-dkms

# 4. Unload conflicting modules
sudo modprobe -r b44 b43 b43legacy ssb brcmsmac bcma

# 5. Load the wl driver
sudo modprobe wl

# 6. Install WiFi tools
sudo apt install wireless-tools wpasupplicant network-manager

# 7. Reboot
sudo reboot
```

**After reboot:**
```bash
# Check if WiFi is working
ip link show

# If WiFi interface shows as "dormant" (common issue), try:
# Put MacBook to sleep and wake it up (close lid, wait 5 sec, open)
# This often activates the WiFi interface

# Configure WiFi
sudo nmtui
```

#### Option 2: Open-Source `b43` Driver (Alternative)

```bash
# 1. Add non-free repositories (same as above)
sudo nano /etc/apt/sources.list
# Add "non-free-firmware" to each line

# 2. Update package list
sudo apt update

# 3. Install b43 firmware installer
sudo apt install firmware-b43-installer

# 4. Load the driver
sudo modprobe b43

# 5. Reboot
sudo reboot
```

**Note**: The `b43` driver may not support 5 GHz networks and can have stability issues. The `wl` driver (Option 1) is recommended for better performance.

**Troubleshooting WiFi "Dormant" State:**
If `ip link show` shows your WiFi interface as "dormant":
```bash
# Try putting MacBook to sleep and waking it
# Close lid, wait 5 seconds, open lid
# WiFi should activate

# Or manually:
sudo ip link set wlp2s0 up  # Replace wlp2s0 with your interface name
```

### Trackpad Configuration

**MacBook trackpads may need configuration:**

```bash
# Install synaptics driver (for older trackpads)
sudo apt install xserver-xorg-input-synaptics

# Or libinput (for newer trackpads, usually pre-installed)
sudo apt install xserver-xorg-input-libinput

# Configure tap-to-click and scrolling
# Edit: /etc/X11/xorg.conf.d/40-libinput.conf (if using libinput)
```

### Graphics Configuration (NVIDIA GeForce 330M)

**Your MacBook Pro 2010 has dual graphics:**
- Intel HD Graphics (integrated) - works out of the box
- NVIDIA GeForce 330M (discrete) - optional, requires drivers

**For server use (website hosting), Intel graphics is sufficient:**

```bash
# Check which graphics are active
lspci | grep -i vga

# Intel graphics should work automatically
# No additional configuration needed for server use
```

**If you want to use NVIDIA graphics (optional, not needed for server):**
```bash
# Install NVIDIA drivers (if needed)
sudo apt install nvidia-driver

# Note: For a server running just a website, this is unnecessary
# Intel graphics uses less power and is perfectly fine
```

### Display Brightness (if needed)

```bash
# Install brightness control
sudo apt install brightnessctl

# Or use xbacklight
sudo apt install xbacklight

# Set brightness (0-100)
sudo brightnessctl set 50%

# For MacBook Pro 2010, you may need to adjust in BIOS/EFI
# or use function keys if desktop environment is installed
```

### Audio Configuration (NVIDIA MCP89)

**Your MacBook Pro 2010 audio should work out of the box:**

```bash
# Install audio utilities
sudo apt install alsa-utils pulseaudio

# Test audio
alsamixer

# If speakers are quiet by default:
# In alsamixer, adjust "Surround" speakers
# Use arrow keys to navigate, 'M' to unmute
```

**Audio should work automatically with the `snd-hda-intel` driver (pre-installed).**

## Website Setup (Next.js)

### Step 1: Install Node.js 20+

```bash
# Add NodeSource repository
curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -

# Install Node.js
sudo apt install -y nodejs

# Verify installation
node --version
npm --version
```

### Step 2: Install Your Website

```bash
# Navigate to your website directory
cd /path/to/website

# Install dependencies
npm install

# Build the website
npm run build

# Test the build
npm start
```

### Step 3: Set Up Production Server

**Option A: Using PM2 (Recommended)**

```bash
# Install PM2 globally
sudo npm install -g pm2

# Start your website
pm2 start npm --name "syste26-website" -- start

# Make PM2 start on boot
pm2 startup
pm2 save
```

**Option B: Using systemd service**

```bash
# Create service file
sudo nano /etc/systemd/system/syste26-website.service
```

Add this content:
```ini
[Unit]
Description=SYSTE-26 Next.js Website
After=network.target

[Service]
Type=simple
User=yourusername
WorkingDirectory=/path/to/website
Environment=NODE_ENV=production
ExecStart=/usr/bin/npm start
Restart=always

[Install]
WantedBy=multi-user.target
```

Then enable and start:
```bash
sudo systemctl daemon-reload
sudo systemctl enable syste26-website
sudo systemctl start syste26-website
```

### Step 4: Set Up Reverse Proxy (Optional but Recommended)

**Install nginx:**

```bash
sudo apt install nginx

# Configure nginx
sudo nano /etc/nginx/sites-available/syste26
```

Add configuration:
```nginx
server {
    listen 80;
    server_name your-domain.com;  # or localhost

    location / {
        proxy_pass http://localhost:3000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_cache_bypass $http_upgrade;
    }
}
```

Enable site:
```bash
sudo ln -s /etc/nginx/sites-available/syste26 /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```

## Energy Optimization

### Install Power Management Tools

```bash
# Install TLP (power management)
sudo apt install tlp tlp-rdw

# Start TLP
sudo tlp start

# Enable TLP on boot
sudo systemctl enable tlp
```

### Configure CPU Governor

```bash
# Install cpufrequtils
sudo apt install cpufrequtils

# Set to powersave mode (for maximum efficiency)
sudo cpufreq-set -g powersave

# Or ondemand (balanced)
sudo cpufreq-set -g ondemand

# Make it permanent
sudo nano /etc/default/cpufrequtils
# Add: GOVERNOR="powersave"
```

### Disable Unnecessary Services

```bash
# List running services
systemctl list-units --type=service --state=running

# Disable unnecessary services (examples)
sudo systemctl disable bluetooth  # If not using Bluetooth
sudo systemctl disable avahi-daemon  # If not using network discovery
```

## Troubleshooting

### WiFi Not Working (MacBook Pro 2010 Specific)

```bash
# Check WiFi card (should show BCM4322)
lspci | grep -i network

# Check if driver is loaded
lsmod | grep -E "wl|b43"

# Check WiFi interface status
ip link show
# If shows "dormant", try sleep/wake workaround

# Install proprietary driver (if not already installed)
sudo apt install broadcom-sta-dkms

# Unload conflicting modules
sudo modprobe -r b44 b43 b43legacy ssb brcmsmac bcma

# Load wl driver
sudo modprobe wl

# If interface is dormant, try:
# 1. Put MacBook to sleep (close lid)
# 2. Wait 5 seconds
# 3. Wake it up (open lid)
# This often activates WiFi

# Or manually bring interface up
sudo ip link set wlp2s0 up  # Replace with your interface name

# Restart networking
sudo systemctl restart NetworkManager
```

### Trackpad Not Working

```bash
# Check if trackpad is detected
xinput list

# Install drivers
sudo apt install xserver-xorg-input-libinput
sudo apt install xserver-xorg-input-synaptics
```

### Display Issues

```bash
# Check graphics
lspci | grep -i vga

# Install drivers (if needed)
# For Intel: Usually works out of the box
# For NVIDIA: sudo apt install nvidia-driver
```

### Boot Issues

```bash
# If GRUB doesn't show up
sudo update-grub

# If can't boot from USB
# Try: Hold Option key longer
# Or: Use rEFInd boot manager
```

## Useful Commands

```bash
# Check system info
uname -a
lscpu
free -h
df -h

# Check network
ip addr
ping google.com

# Check services
systemctl status nginx
systemctl status syste26-website

# View logs
journalctl -u syste26-website -f
tail -f /var/log/nginx/error.log
```

## Next Steps

1. ✅ Install Debian 12
2. ✅ Set up WiFi
3. ✅ Install Node.js
4. ✅ Deploy your website
5. ✅ Configure power management
6. ✅ Set up reverse proxy (optional)
7. ✅ Configure firewall (optional but recommended)

## Security Recommendations

```bash
# Install and configure firewall
sudo apt install ufw
sudo ufw allow 22/tcp  # SSH
sudo ufw allow 80/tcp  # HTTP
sudo ufw allow 443/tcp  # HTTPS
sudo ufw enable

# Disable root SSH login (edit /etc/ssh/sshd_config)
# Set: PermitRootLogin no
sudo systemctl restart sshd
```

## MacBook Pro 2010 17-inch Specific Resources

- **Debian Wiki for MacBook Pro 7,1** (your model): https://wiki.debian.org/InstallingDebianOn/Apple/MacBookPro/7-1
- **Broadcom WiFi Drivers**: https://wiki.debian.org/bcm43xx
- **Debian Installation Guide**: https://www.debian.org/releases/stable/installmanual
- **Node.js Installation**: https://github.com/nodesource/distributions
- **PM2 Documentation**: https://pm2.keymetrics.io/

## Quick Reference for Your MacBook Pro 2010

**Hardware Summary:**
- ✅ **CPU**: Intel Core i7 (64-bit) - Use amd64 ISO
- ⚠️ **WiFi**: Broadcom BCM4322 - Install `broadcom-sta-dkms` driver
- ✅ **Graphics**: Intel HD (works) + NVIDIA 330M (optional)
- ✅ **Audio**: NVIDIA MCP89 - Works with snd-hda-intel
- ✅ **Bluetooth**: Broadcom BCM2046 - Works with btusb
- ✅ **Suspend**: Works out of the box

**Installation Priority:**
1. Use firmware ISO (includes WiFi drivers)
2. Install `broadcom-sta-dkms` for WiFi after installation
3. If WiFi is "dormant", use sleep/wake workaround
4. For server use, Intel graphics is sufficient (no NVIDIA needed)

---

**Good luck with your installation!** 🚀

Your 2010 MacBook Pro 17-inch should work well with Debian 12. The main thing to watch for is the WiFi driver installation - make sure to install `broadcom-sta-dkms` for best results.

